function writeQ(DataStruct, varargin)
% writeQ - Saves any struct from IsoNQ as CSV, XLSX, MAT, or TXT files in a named subdirectory.
%
% Inputs:
%   DataStruct - Struct containing data to be saved. [struct]
%
% Optional Parameters:
%   'OutputDir'  - Parent directory where files will be saved. [char/string]
%                  Default: './output'
%   'Format'     - Output file format: 'csv', 'xlsx', 'mat', or 'txt'. Default: 'xlsx'
%   'Timestamp'  - Whether to include a timestamped directory (true/false). Default: true
%
% Example:
%   writeQ(PX, 'OutputDir', './output', 'Format', 'xlsx', 'Timestamp', false);
%   % Saves files in './output/PX/'
%
%   writeQ(QX, 'OutputDir', './results', 'Format', 'mat', 'Timestamp', true);
%   % Saves QX as a .mat file inside './results/QX_YYYYMMDD_HHMMSS/'

% --- Define Expected Inputs and Defaults ---
validStruct = @(x) isstruct(x);
validDir = @(x) ischar(x) || isstring(x);
validFormat = @(x) any(validatestring(x, {'csv', 'xlsx', 'mat', 'txt'}));
validTimestamp = @(x) islogical(x) || isnumeric(x);

% Default values
defaultFormat = 'xlsx';
defaultOutputDir = 'output';
defaultTimestamp = true;

% Set up input parser
p = inputParser;
addRequired(p, 'DataStruct', validStruct);
addParameter(p, 'OutputDir', defaultOutputDir, validDir);
addParameter(p, 'Format', defaultFormat, validFormat);
addParameter(p, 'Timestamp', defaultTimestamp, validTimestamp);

% Parse inputs
parse(p, DataStruct, varargin{:});

% Extract parsed values
OutputDir = p.Results.OutputDir;
Format = p.Results.Format;
UseTimestamp = p.Results.Timestamp;

% --- Ensure the Output Directory is Named After the Struct ---
structName = inputname(1); % Extract variable name (e.g., 'PX')

if isempty(structName)
    structName = 'IsoNQ_Output'; % Default if struct is anonymous
end

% Generate timestamped directory if enabled
if UseTimestamp
    timestamp = datestr(now, 'yyyymmdd_HHMMSS'); % Format: YYYYMMDD_HHMMSS
    fullOutputDir = fullfile(OutputDir, sprintf('%s_%s', structName, timestamp));
else
    fullOutputDir = fullfile(OutputDir, structName);
end

% Create the directory if it does not exist
if ~exist(fullOutputDir, 'dir')
    mkdir(fullOutputDir);
end

% Initialize logging
logFileName = fullfile(fullOutputDir, 'writeQ_log.txt');
logFileID = fopen(logFileName, 'w');
fprintf(logFileID, "IsoNQ writeQ Log - %s\n\n", datestr(now));

% If saving as .mat, store the entire struct in one file
if strcmp(Format, 'mat')
    fileName = fullfile(fullOutputDir, sprintf('%s.mat', structName));
    save(fileName, '-struct', 'DataStruct');
    fprintf(logFileID, "Saved: %s\n", fileName);
    fprintf('Struct saved as %s\n', fileName);
    fclose(logFileID);
    return;
end

% --- Extract Fields from Struct and Save Data ---
fieldNames = fieldnames(DataStruct);
xlsxFileName = fullfile(fullOutputDir, sprintf('%s.xlsx', structName));

for i = 1:length(fieldNames)
    fieldName = fieldNames{i};
    data = DataStruct.(fieldName);

    % Construct file path (except for xlsx)
    fileName = fullfile(fullOutputDir, sprintf('%s.%s', fieldName, Format));

    % Handle different data types
    if strcmp(Format, 'xlsx')
        try
            if istable(data)
                % Detect and unnest nested tables before writing
                nestedVars = varfun(@(x) istable(x), data, 'OutputFormat', 'uniform');
                if any(nestedVars)
                    data = splitvars(data);
                end
                writetable(data, xlsxFileName, 'Sheet', fieldName);
            else
                writematrix(data, xlsxFileName, 'Sheet', fieldName);
            end
            fprintf(logFileID, "Saved to XLSX (Sheet: %s)\n", fieldName);
        catch ME
            warning('Could not write field "%s" due to error: %s', fieldName, ME.message);
        end
    elseif istable(data)
        try
            if strcmp(Format, 'csv')
                writetable(data, fileName);
            else
                writetable(data, fileName, 'FileType', 'spreadsheet');
            end
        catch ME
            warning('Skipping field "%s" - Unable to write due to error: %s', fieldName, ME.message);
        end
    elseif isnumeric(data) || iscell(data)
        tableData = array2table(data);
        if strcmp(Format, 'csv')
            writetable(tableData, fileName);
        else
            writetable(tableData, fileName, 'FileType', 'spreadsheet');
        end
    elseif ischar(data) || isstring(data)
        fileID = fopen(fileName, 'w');
        fprintf(fileID, '%s\n', data);
        fclose(fileID);
    elseif strcmp(Format, 'txt')
        fileID = fopen(fileName, 'w');
        fprintf(fileID, 'Field Name: %s\n\n', fieldName);
        disp(data, fileID);
        fclose(fileID);
    else
        warning('Skipping field "%s" - Unsupported data type.', fieldName);
    end

    fprintf(logFileID, "Saved: %s\n", fileName);
end

fclose(logFileID);
fprintf('Files saved in %s\n', fullOutputDir);
end
