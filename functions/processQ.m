function [PX] = processQ(fnames,IX,varargin)
% processQ - Executes the full IsoNQ data processing workflow.
%
% This function performs parsing, outlier removal, drift correction, size correction, 
% scale normalization, uncertainty propagation, and weighted averaging to generate final 
% processed isotope data. Optionally, it can generate plots and log processing parameters.
%
% Inputs:
%   fnames      - List of file names containing raw Qtegra CSV data. [Fx1 string array]
%                 Each file should be exported using Qtegra’s CSV format without modifications.
%   IX          - Structured input containing all processing parameters, as created by infoQ. [struct]
%   Plot        - Option to generate and save figures to the output directory ('yes' or 'no'). [char]
%                 Default: 'no'
%
% Outputs:
%   PX          - Struct containing all processed data, metadata, and results. [struct]
%                 Includes raw and corrected isotope values, propagated uncertainties, 
%                 and final weighted means. See PX.Descriptions for details on individual fields.
%
% Notes:
% - This function applies all data corrections sequentially, including:
%   (1) Parsing input files and extracting sample/reference gas peaks.
%   (2) Performing outlier analysis **(set in infoQ; options: 'mad', 'zscore', 'iqr', 'grubbs', 'none')**.
%   (3) Applying compound-specific or global drift correction **(set in infoQ; options: 'CompoundSpecific', 'GlobalDrift')**.
%   (4) Applying size correction.
%   (5) Applying regression-based or two-point scale correction **(set in infoQ; options: 'Regression', 'TwoPoint')**.
%   (6) Calculating uncertainty propagation from sample and reference gas.
%   (7) Generating final weighted means for selected compounds.
% - All intermediate outputs are stored in PX and can be accessed individually if needed.
% - If 'Plot' is set to 'yes', figures will be generated and saved into an output directory called 'fig'.
%
% Example Usage:
%   IX = infoQ(...);  % Define input parameters
%   fnames = ["sample1.csv", "sample2.csv"];
%   PX = processQ(fnames, IX, 'Plot', 'yes');
%
% See Also: infoQ, parseQ, correctQ, errorQ, weightQ

tic

defPlot = 'no';

expPlot = {'yes','no'};

p = inputParser;
validfnames = @(x) isstring(x);
validIX = @(x) isstruct(x);
validPlot = @(x) any(validatestring(x,expPlot));

addRequired(p,'fnames',validfnames)
addRequired(p,'IX',validIX)
addParameter(p,'Plot',defPlot,validPlot)

parse(p,fnames,IX,varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

fnames = p.Results.fnames;
IX = p.Results.IX;
Plot = char(p.Results.Plot);

Mode = IX.Mode;
RemoveFlags = IX.RemoveFlags;

[EX] = parseQ(fnames,'Mode',Mode);
[AX] = findQ(EX,IX);
[QX] = correctQ(AX,IX,'Plot',Plot);
[CX] = combineQ(QX,'Plot',Plot);
[AQ] = assignQ(CX,IX,'Plot',Plot);

[UCX] = errorQ(AQ,AX,IX,'Plot',Plot);
[WT] = weightQ(AQ,IX,'Error',UCX,'Plot',Plot,'RemoveFlags',RemoveFlags);

PX = struct();
PX.MeanDelta = AQ.MeanDelta;
PX.SDDelta = AQ.SDDelta;
PX.Weighted = WT;
PX.TotalUncertainty = UCX.TotalUncertainty;
PX.RefgasUncertainty = UCX.RefgasUncertainty;
PX.SampleUncertainty = UCX.SampleUncertainty;
PX.FlaggedUncertainty = UCX.FlaggedUncertainty;
PX.FlaggedRemovedDelta = UCX.FlaggedRemovedDelta;
PX.FlaggedRemovedUncertainty = UCX.FlaggedRemovedUncertainty;
PX.FlaggedSamples = UCX.FlaggedSamples;
PX.Covariance = UCX.Covariance;
PX.Standards = CX.Standard;
PX.Samples = CX.Sample;
PX.FileName = fnames;

vn = {'Field','Description'};
Field = string(fieldnames(PX));
Description = string({"Mean isotope values of measured samples for each compound.",...
    "Standard deviation of isotope values for measured samples for each compound.",...
    "Weighted mean isotope values across replicates.",...
    "Total analytical uncertainty, combining reference gas and sample contributions.",...
    "Analytical uncertainty associated with the reference gas.",...
    "Analytical uncertainty associated with the sample.",...
    "Table of total uncertainty values exceeding the flag threshold; values below threshold are removed.",...
    "Table of isotope values with flagged samples removed.",...
    "Table of total uncertainty values with flagged samples removed.",...
    "List of flagged samples, including the number of flagged compounds and chain lengths.",...
    "Covariance of reference gas and sample uncertainties.",...
    "Table of isotope values and metadata for lab standards.",...
    "Table of isotope values and metadata for samples.",...
    "File names of processed Qtegra data files."})';
DT = array2table(zeros(length(Field),2));
DT.Properties.VariableNames = vn;
DT.Field = Field;
DT.Description = Description;
PX.Descriptions = DT;

[PX.Function] = deal('processQ');

toc
end