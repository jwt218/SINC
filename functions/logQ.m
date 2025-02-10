function logQ(IX, OutputDir)

% Create timestamped filename
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
logFile = fullfile(OutputDir, 'log', ['IsoNQ_Parameters_' timestamp '.txt']);

% Ensure log directory exists
if ~exist(fullfile(OutputDir, 'log'), 'dir')
    mkdir(fullfile(OutputDir, 'log'));
end

% Open file for writing
fid = fopen(logFile, 'w');

% Write header
fprintf(fid, 'IsoNQ Parameter Log - %s\n\n', datestr(now));

% Write each parameter from IX struct
fields = fieldnames(IX);
for i = 1:length(fields)
    value = IX.(fields{i});

    % Convert numeric arrays to string
    if isnumeric(value)
        valueStr = mat2str(value);
    else
        valueStr = string(value);
    end

    fprintf(fid, '%s: %s\n', fields{i}, valueStr);
end

% Close file
fclose(fid);
%fprintf('Parameter log saved: %s\n', logFile);
end
