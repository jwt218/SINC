function MX = mergeQ(varargin)
% mergeQ - Merges multiple IsoNQ output structs (PX) into a combined dataset.
%
% Inputs:
%   Multiple PX structs from processQ (e.g., PX1, PX2, PX3)
%
% Outputs:
%   MX - A merged struct containing combined sample and standard data, filling missing values.
%
% Example:
%   MX = mergeQ(PX1, PX2, PX3);

% Ensure at least two inputs
if nargin < 2
    error('mergeQ requires at least two PX structs as inputs.');
end

% Initialize MX with the first PX struct
MX = varargin{1};

% Get list of expected fields in PX
expectedFields = fieldnames(MX);

% Loop through additional datasets and merge each field
for i = 2:nargin
    PX = varargin{i};

    % Loop through each field in PX
    for j = 1:length(expectedFields)
        fieldName = expectedFields{j};

        if isfield(PX, fieldName)
            if isfield(MX, fieldName)
                % If both are tables, perform ID-based merge and fill missing values
                if istable(MX.(fieldName)) && istable(PX.(fieldName))
                    joinKey = intersect(MX.(fieldName).Properties.VariableNames, PX.(fieldName).Properties.VariableNames);
                    
                    if any(strcmp(joinKey, "ID")) || any(strcmp(joinKey, "Identifier"))
                        keyColumn = joinKey{1}; % Use the matching key
                        mergedTable = outerjoin(MX.(fieldName), PX.(fieldName), ...
                            'Keys', keyColumn, 'MergeKeys', true, 'Type', 'left');

                        % Convert to array for isnan check
                        numericData = table2array(mergedTable(:, 2:end));

                        % Find missing values
                        missingMask = isnan(numericData);

                        % Convert PX data to numeric for replacement
                        pxDataArray = table2array(PX.(fieldName)(:, 2:end));

                        % Fill missing values
                        numericData(missingMask) = pxDataArray(missingMask);

                        % Convert back to table and store in MX
                        mergedTable(:, 2:end) = array2table(numericData);
                        MX.(fieldName) = mergedTable;
                    end
                elseif isstruct(MX.(fieldName)) && isstruct(PX.(fieldName))
                    % Handle nested subfields (e.g., PX.Weighted.Data)
                    subfields = fieldnames(MX.(fieldName));
                    for k = 1:length(subfields)
                        subfieldName = subfields{k};
                        if isfield(PX.(fieldName), subfieldName)
                            MX.(fieldName).(subfieldName) = mergeQ(MX.(fieldName).(subfieldName), PX.(fieldName).(subfieldName));
                        end
                    end
                elseif isnumeric(MX.(fieldName)) && isnumeric(PX.(fieldName))
                    % Fill missing values in numeric arrays
                    missingIdx = isnan(MX.(fieldName));
                    MX.(fieldName)(missingIdx) = PX.(fieldName)(missingIdx);
                end
            else
                % If field doesn't exist in MX, add it
                MX.(fieldName) = PX.(fieldName);
            end
        end
    end
end

% Merge FileNames correctly
if ischar(MX.FileName) || isstring(MX.FileName)
    MX.FileName = {MX.FileName};
end
MX.FileName = unique([MX.FileName; PX.FileName]);

% Maintain function metadata
[MX.Function] = deal('mergeQ');

fprintf('Successfully merged %d IsoNQ datasets with missing-value filling and nested struct handling.\n', nargin);
end
