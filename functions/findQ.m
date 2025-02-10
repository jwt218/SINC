function [AX] = findQ(EX, IX, varargin)
% findQ - Distributes analysis (sample and standard) and reference gas data
%         into separate tables based on SampleRT and RefgasRT.
%
% Inputs:
%   EX       - Output from parseQ [struct]
%   IX       - Output from infoQ [struct]
%              Required input parameters from IX include SampleRT and
%              RefgasRT. Defaults will be used for Compounds, AdditionalComp, KnownValues,
%              KnownSD, PeakWindow, RefgasWindow, and ReferenceShift if not specified.
%
% Outputs:
%   AX       - Struct with fields:
%              AX.Analysis - Processed sample and standard data
%              AX.Refgas   - Extracted reference gas data
%              AX.FileName - Associated file names
%
% Example Usage:
%   AX = findQ(EX, IX);

%% --- Input Parsing and Validation ---
p = inputParser;
validEX = @(x) isstruct(x);
validIX = @(x) isstruct(x);

addRequired(p, 'EX', validEX)
addRequired(p, 'IX', validIX)

parse(p, EX, IX, varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

EX = p.Results.EX(:);
IX = p.Results.IX(:);

%% --- Extract Necessary Parameters from IX ---
HX = EX;
sampleRT = IX.SampleRT;
refgasRT = IX.RefgasRT;
compounds = [IX.Compounds; IX.AdditionalComp];  % Combined list of compounds
knownValues = IX.KnownValues;
knownSD = IX.KnownSD;
peakWindow = IX.PeakWindow;
refgasWindow = IX.RefgasWindow;
referenceShift = IX.ReferenceShift;

AX = struct();

%% --- Process Each Data File ---
for j = 1:length(EX)
    H = HX(j).Data;
    Ert = H.RetentionTime;
    vns = H.Properties.VariableNames;
    vts = {'string', 'double', 'double', 'double', 'double', ...
           'double', 'double', 'double'};

    sz = size(H);
    ll = sz(1);
    A = table('Size', sz, 'VariableTypes', vts);
    A.Properties.VariableNames = vns;
    cnum = length(A.Properties.VariableNames);
    A(:, 2:cnum) = array2table(repmat(9999, ll, cnum-1));
    A.KnownValue = repmat(9999, ll, 1);
    A.KnownSD = repmat(9999, ll, 1);
    A.Component = repmat(9999, ll, 1);

    %% --- Identify and Extract Sample Peaks ---
    dif = peakWindow;
    for i = 1:length(sampleRT)
        if ~isnan(sampleRT(i))
            ki = find(Ert >= sampleRT(i) - dif & Ert <= sampleRT(i) + dif);
            A(ki, 1:cnum) = H(ki, 1:cnum);
            if i <= length(knownValues)
                A.KnownValue(ki) = knownValues(i);
                A.KnownSD(ki) = knownSD(i);
                A.Component(ki) = compounds(i);
            else
                A.KnownValue(ki) = NaN;
                A.KnownSD(ki) = NaN;
                A.Component(ki) = compounds(i);
            end
        end
    end

    % Remove unfilled rows
    A(A.(3) == 9999, :) = [];

    % Compute cumulative retention time
    A.CumulativeRT = cumsum(table2array(A(:, "RetentionTime")));

    %% --- Identify and Extract Reference Gas Peaks ---
    Ag = table('Size', sz, 'VariableTypes', vts);
    Ag.Properties.VariableNames = vns;
    Ag(:, 2:cnum) = array2table(repmat(9999, ll, cnum-1));

    dif = refgasWindow;
    for i = 1:length(refgasRT)
        if ~isnan(refgasRT(i))
            ki = find(Ert >= refgasRT(i) - dif & Ert <= refgasRT(i) + dif);
            Ag(ki, 1:cnum) = H(ki, 1:cnum);
            if i <= length(knownValues)
                Ag.KnownValue(ki) = knownValues(i);
                Ag.Component(ki) = compounds(i);
            else
                Ag.KnownValue(ki) = NaN;
                Ag.Component(ki) = compounds(i);
            end
        end
    end

    % Remove unfilled rows
    Ag(Ag.(3) == 9999, :) = [];

    % Apply reference gas shift correction
    A.RawDeltaValue = A.RawDeltaValue - referenceShift;
    Ag.RawDeltaValue = Ag.RawDeltaValue - referenceShift;

    %% --- Store Processed Data ---
    AX(j).Analysis = A;
    AX(j).Refgas = Ag;
    AX(j).FileName = HX(j).FileName;
end

%% --- Assign Function Name for Reference ---
[AX.Function] = deal('findQ');

end
