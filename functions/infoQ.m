function [IX] = infoQ(varargin)
% infoQ - Organizes input parameters into a structured format for IsoNQ processing.
%
% This function initializes and validates user-defined parameters for the IsoNQ processing pipeline.
% Parameters include retention times, compound selection, drift/scale correction methods, uncertainty thresholds,
% and outlier detection settings.
%
% Inputs (Optional Parameters with Defaults):
%   'SampleRT'          - Retention times for sample peaks (in measurement units). [Nx1 double], Default: []
%   'RefgasRT'          - Retention times for reference gas injections. [Mx1 double], Default: []
%   'Compounds'         - Primary list of compound numbers to analyze (e.g., [29 31] for C29 and C31). [Px1 double], Default: 16:30
%   'AdditionalComp'    - Additional compounds to include in analysis. [Qx1 double], Default: []
%   'KnownValues'  - Known isotope values for reference compounds.
%                    - Can be a **vector** of known values (e.g., [-26.1, -31.8, -32.7]).
%                    - Can be a **standard name** (e.g., 'A6', 'B4') to apply built-in values.
%                    - Can be **'match'**, which automatically follows StandardID.
%                    - If a **standard name** is selected, values are chosen based on 'Mode' (δ¹³C for 'C', δD for 'H').
%                    Default: 'match'
%   'KnownSD'      - Standard deviations of known isotope values.
%                    - If 'KnownValues' is a **vector**, this must also be manually set as a vector.
%                    - If 'KnownValues' is a **standard name** or 'match', KnownSD is assigned automatically.
%                    Default: 'match'
%   'PeakWindow'        - Detection window for sample peaks (measurement units). [1x1 double], Default: 5
%   'RefgasWindow'      - Detection window for identifying reference gas peaks. [1x1 double], Default: 15
%   'SampleMatchLen'    - Number of leading characters used to match sample replicates (e.g., 3 for IDs like SM1_a and SM1_b). [1x1 double], Default: 3
%   'AgeFile'           - File path for sample age data, if available. [char], Default: []
%   'WeightedComp'      - List of compound numbers to use for weighted averaging. [Wx1 double], Default: [25 27 29]
%   'UncertaintyThreshold' - Threshold for removing high-uncertainty samples when 'RemoveFlags' is enabled. [1x1 double], Default: 4
%   'PolyOrder'         - Order of polynomial used for drift correction (1 = linear, 2 = quadratic). [1x1 double], Default: 2
%   'RemoveFlags'       - Whether to remove flagged samples exceeding uncertainty threshold ('yes' or 'no'). [char], Default: 'no'
%   'Mode'              - Isotope analysis mode: 'C' for carbon, 'H' for hydrogen. [char], Default: 'C'
%   'ReferenceShift'    - Uniform shift applied to all values to adjust reference gas offset. [1x1 double], Default: 0
%   'DriftCorrection'   - Method used for drift correction: 'CS' (compound-specific) or 'Bulk'. [char], Default: 'CS'
%   'ScaleCorrection'   - Scale correction method: 'Global' (all compounds) or '2PT' (two-point correction). [char], Default: 'Global'
%   'Scale2PTComp'      - Compounds used for two-point scale correction. [2x1 double], Default: [23 29]
%   'OutlierMethod'     - Method for outlier detection before applying corrections. Options:
%                         'mad'    - Median Absolute Deviation (Default, robust)
%                         'zscore' - Z-score (removes values >3 standard deviations)
%                         'iqr'    - Interquartile Range (removes values outside 1.5×IQR range)
%                         'grubbs' - Grubbs' test (detects single extreme outliers in small datasets)
%                         'none'   - No outlier removal
%   'OutlierThreshold'  - Multiplier for thresholding outliers (applicable to MAD and Z-score methods). [1x1 double], Default: 5
%
% Example Usage:
%   IX = infoQ('SampleRT', [1 3 5 6], 'Compounds', [25 26 27 28], 'StandardID', 'A6', 'Mode', 'C', ...
%              'OutlierMethod', 'iqr', 'OutlierThreshold', 4);
%
% Outputs:
%   IX - Struct containing all validated input parameters.
%
% Notes:
%   - The default **outlier detection method** is 'mad' (Median Absolute Deviation).
%   - OutlierThreshold is only applied to the **'mad'** and **'zscore'** methods.
%   - To disable outlier removal, set 'OutlierMethod' to 'none'.
%   - Outlier detection is performed **before** applying drift, size, and scale corrections.
%   - If 'KnownValues' is **'match'**, it automatically follows StandardID.
%   - If 'KnownValues' is a **standard name**, KnownSD is also automatically assigned.
%   - If 'KnownValues' is a **custom vector**, KnownSD must be manually provided.
%   - Standard names available: 'A6', 'B4' (more can be added).
%   - The **Mode ('C' or 'H')** determines whether **δ¹³C or δD** values are used.



%% --- Define Default Values ---

defSampleRT = [];
defRefgasRT = [];
defCompounds = 16:30;
defAdditionalComp = [];
defStandardID = "A6";
defPeakWindow = 5;
defRefgasWindow = 15;
defSampleMatchLen = 3;
defAgeFile = [];
defWeightedComp = [25 27 29];
defUncertaintyThreshold = 4;
defPolyOrder = 2;
defRemoveFlags = 'no';
defMode = 'C';
defReferenceShift = 0;
defDriftCorrection = 'CS';
defScaleCorrection = 'Global';
defScale2PTComp = [23 29];
defOutlierMethod = 'mad';
defOutlierThreshold = 5;
defLogParams = 'yes';
defLogDir = './output';

% --- Predefined Known Values and SD for Standards ---
predefinedStandards = struct();

% δ¹³C Values
predefinedStandards.A6.C.Values = [-26.15 -31.88 -32.7 -31.99 -33.97 -28.83 -33.77 -33.37 -32.13 -28.46 -32.94 -30.49 -33.2 -29.1 -29.84];
predefinedStandards.A6.C.SD = [0.02 0.02 0.01 0.01 0.02 0.02 0.02 0.03 0.02 0.02 0.01 0.01 0.01 0.01 0.01];

% δD Values
predefinedStandards.A6.H.Values = [-150.2 -145.8 -160.3 -155.9 -162.7 -148.1 -158.9 -153.4 -161.2 -147.6 -154.9 -157.3 -163.5 -151.4 -159.8];
predefinedStandards.A6.H.SD = [2.0 2.1 1.8 2.0 2.2 2.3 1.9 2.4 2.1 2.0 1.7 2.0 1.9 2.1 1.8];

% δ¹³C for B4
predefinedStandards.B4.C.Values = [-30.10 -33.05 -34.20 -32.80 -31.60 -28.70 -32.90 -30.50 -33.00 -29.80 -30.40 -32.00 -34.50 -29.90 -31.20];
predefinedStandards.B4.C.SD = [0.02 0.02 0.01 0.01 0.02 0.02 0.02 0.03 0.02 0.02 0.01 0.01 0.01 0.01 0.01];

% δD for B4
predefinedStandards.B4.H.Values = [-145.2 -143.0 -159.2 -154.7 -160.5 -149.0 -155.5 -150.3 -162.1 -148.3 -151.9 -155.6 -164.0 -152.5 -157.9];
predefinedStandards.B4.H.SD = [2.1 2.0 1.9 2.1 2.3 2.2 1.8 2.5 2.0 1.9 1.7 2.2 1.9 2.0 1.8];

defKnownValues = 'match';  % Default to auto-match StandardID
defKnownSD = 'match';      % Default to auto-match StandardID SD

expRemoveFlags = {'yes','no'};
expMode = {'C','H'};
expDriftCorrMethod = {'CS','Bulk'};
expScaleCorrMethod = {'Global','Scale2PT'};
expOutlierMethod = {'mad', 'zscore', 'iqr', 'grubbs', 'none'};
expLogParams = {'yes', 'no'};

p = inputParser;
validSampleRT = @(x) isa(x,'double') && isvector(x);
validRefgasRT = @(x) isa(x,'double') && isvector(x);
validCompounds = @(x) isa(x,'double') && isvector(x);
validAdditionalComp = @(x) isa(x,'double');
validStandardID = @(x) isstring(x) || ischar(x);
validPeakWindow = @(x) isa(x,'double') && isscalar(x);
validRefgasWindow = @(x) isa(x,'double') && isscalar(x);
validSampleMatchLen = @(x) isa(x,'double') && isscalar(x);
validAgeFile = @(x) isstring(x) || ischar(x);
validWeightedComp = @(x) isa(x,'double') && length(x) > 1;
validUncertaintyThreshold = @(x) isa(x,'double') && isscalar(x);
validPolyOrder = @(x) isa(x,'double') && isscalar(x);
validRemoveFlags = @(x) any(validatestring(x,expRemoveFlags));
validMode = @(x) any(validatestring(x,expMode));
validReferenceShift = @(x) isa(x,'double') && isscalar(x);
validDriftCorrection = @(x) any(validatestring(x,expDriftCorrMethod));
validScaleCorrection = @(x) any(validatestring(x,expScaleCorrMethod));
validScale2PTComp = @(x) isa(x,'double') && length(x) == 2;
validOutlierMethod = @(x) any(validatestring(x,expOutlierMethod));
validOutlierThreshold = @(x) isa(x,'double') && isscalar(x);
validLogParams =  @(x) any(validatestring(x, expLogParams));
validLogDir = @(x) isstring(x) || ischar(x);

addParameter(p,'SampleRT',defSampleRT,validSampleRT)
addParameter(p,'RefgasRT',defRefgasRT,validRefgasRT)
addParameter(p,'Compounds',defCompounds,validCompounds)
addParameter(p,'AdditionalComp',defAdditionalComp,validAdditionalComp)
addParameter(p,'KnownValues',defKnownValues)
addParameter(p,'KnownSD',defKnownSD)
addParameter(p,'StandardID',defStandardID,validStandardID)
addParameter(p,'PeakWindow',defPeakWindow,validPeakWindow)
addParameter(p,'RefgasWindow',defRefgasWindow,validRefgasWindow)
addParameter(p,'SampleMatchLen',defSampleMatchLen,validSampleMatchLen)
addParameter(p,'AgeFile',defAgeFile,validAgeFile)
addParameter(p,'WeightedComp',defWeightedComp,validWeightedComp)
addParameter(p,'UncertaintyThreshold',defUncertaintyThreshold,validUncertaintyThreshold)
addParameter(p,'PolyOrder',defPolyOrder,validPolyOrder)
addParameter(p,'Mode',defMode,validMode)
addParameter(p,'RemoveFlags',defRemoveFlags,validRemoveFlags)
addParameter(p,'ReferenceShift',defReferenceShift,validReferenceShift)
addParameter(p,'DriftCorrection',defDriftCorrection,validDriftCorrection)
addParameter(p,'ScaleCorrection',defScaleCorrection,validScaleCorrection)
addParameter(p,'Scale2PTComp',defScale2PTComp,validScale2PTComp)
addParameter(p,'OutlierMethod',defOutlierMethod,validOutlierMethod)
addParameter(p,'OutlierThreshold',defOutlierThreshold,validOutlierThreshold)
addParameter(p,'LogParams',defLogParams,validLogParams)
addParameter(p,'LogDir',defLogDir,validLogDir)

%% --- Parse Inputs ---
parse(p,varargin{:})

if ~isempty(fieldnames(p.Unmatched))
    disp('Extra inputs:')
    disp(p.Unmatched)
end

SampleRT = p.Results.SampleRT(:);
RefgasRT = p.Results.RefgasRT(:);
Compounds = p.Results.Compounds(:);
AdditionalComp = p.Results.AdditionalComp(:);
StandardID = p.Results.StandardID;
PeakWindow = p.Results.PeakWindow(:);
RefgasWindow = p.Results.RefgasWindow(:);
SampleMatchLen = p.Results.SampleMatchLen(:);
AgeFile = p.Results.AgeFile;
WeightedComp = p.Results.WeightedComp(:);
UncertaintyThreshold = p.Results.UncertaintyThreshold(:);
PolyOrder = p.Results.PolyOrder(:);
Mode = p.Results.Mode;
RemoveFlags = p.Results.RemoveFlags;
ReferenceShift = p.Results.ReferenceShift(:);
DriftCorrection = p.Results.DriftCorrection;
ScaleCorrection = p.Results.ScaleCorrection;
Scale2PTComp = p.Results.Scale2PTComp;
OutlierMethod = p.Results.OutlierMethod;
OutlierThreshold = p.Results.OutlierThreshold;
LogParams = p.Results.LogParams;
LogDir = p.Results.LogDir;

%% --- Assign Known Values and SD Based on User Input ---
if strcmp(p.Results.KnownValues, 'match')
    stdName = StandardID;
else
    stdName = p.Results.KnownValues;
end

if isfield(predefinedStandards, stdName) && isfield(predefinedStandards.(stdName), Mode)
    KnownValues = predefinedStandards.(stdName).(Mode).Values;
    KnownSD = predefinedStandards.(stdName).(Mode).SD;
else
    error("Standard '%s' does not support Mode '%s'.", stdName, Mode);
end


msg1 = sprintf("'SampleRT' length (%d) must match 'Compounds' + 'AdditionalComp' length (%d).", length(SampleRT), length(Compounds) + length(AdditionalComp));
assert(length(SampleRT) == length(Compounds) + length(AdditionalComp), msg1);

msg2 = sprintf("'Report' length (%d) must match 'ReportSD' length (%d).", length(KnownValues), length(KnownSD));
assert(length(KnownValues) == length(KnownSD), msg2);

%% --- Construct Output Struct ---

IX.SampleRT = SampleRT;
IX.RefgasRT = RefgasRT;
IX.Compounds = Compounds;
IX.AdditionalComp = AdditionalComp;
IX.KnownValues = KnownValues;
IX.KnownSD = KnownSD;
IX.StandardID = string(StandardID);
IX.PeakWindow = PeakWindow;
IX.RefgasWindow = RefgasWindow;
IX.SampleMatchLen = SampleMatchLen;
IX.AgeFile = AgeFile;
IX.WeightedComp = WeightedComp;
IX.UncertaintyThreshold = UncertaintyThreshold;
IX.PolyOrder = PolyOrder;
IX.Mode = Mode;
IX.RemoveFlags = RemoveFlags;
IX.ReferenceShift = ReferenceShift;
IX.DriftCorrection = DriftCorrection;
IX.ScaleCorrection = ScaleCorrection;
IX.Scale2PTComp = Scale2PTComp;
IX.OutlierMethod = OutlierMethod;
IX.OutlierThreshold = OutlierThreshold;
IX.LogParams = LogParams;
IX.LogDir = LogDir;

[IX.Function] = deal('infoQ');

if strcmp(p.Results.LogParams, 'yes')
    logQ(IX, LogDir);
end

end