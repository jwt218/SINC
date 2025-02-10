function [SX] = standardsQ(IX, varargin)
% standardsQ - Extracts and compiles standard data across multiple parsed runs using findQ.
%
% Inputs:
%   IX         - Output from infoQ [struct], used to extract 'StandardID'
%   varargin   - Multiple parsed outputs from parseQ (e.g., EX1, EX2, EX3, ...)
%                Each input must be a struct from parseQ.
%   Plot       - Option to generate diagnostic figures ('yes' or 'no'). Default: 'yes'
%
% Outputs:
%   SX         - Struct with fields:
%                SX.StandardData   - Table of all standard measurements
%                SX.RefgasData     - Table of reference gas measurements
%                SX.FileNames      - List of input files used
%
% Example Usage:
%   SX = standardsQ(IX, EX1, EX2, 'Plot', 'yes');
%
% Notes:
%   - Uses findQ to extract **standard data** from multiple parsed runs.
%   - Generates **diagnostic figures** for drift, residuals, and variability.
%   - The 'StandardID' from 'IX' determines which rows are classified as standards.


%% --- Input Parsing ---
defPlot = 'yes';
expPlot = {'yes', 'no'};

p = inputParser;
addRequired(p, 'IX', @isstruct);  % Ensure IX is a struct
addParameter(p, 'Plot', defPlot, @(x) any(validatestring(x, expPlot)));

% Parse only IX and named parameters first (not varargin yet)
parse(p, IX, varargin{cellfun(@ischar, varargin)});

% Extract named parameters
Plot = char(p.Results.Plot);

% Extract EX structs from varargin
EXList = varargin(~cellfun(@ischar, varargin));  % Get only structs
assert(all(cellfun(@isstruct, EXList)), "All extra inputs must be parsed outputs from parseQ.");

%% --- Extract StandardID from IX ---
StandardID = IX.StandardID;
Mode = IX.Mode;

%% --- Initialize Output Struct ---
SX = struct();
StandardData = table();
RefgasData = table();
FileNames = {};

%% --- Process Each EX Struct with findQ ---
for i = 1:length(EXList)
    EX = EXList{i};
    Q = [];
    for j = 1:length(EX)
        Q = [Q;EX(j).FileName];
    end
    FileNames{i} = Q;

    % Run findQ to extract standard and refgas data
    AX = findQ(EX, IX);
    axy = AX.Analysis;
    axr = AX.Refgas;

    % Append extracted standard and reference gas data
    StandardData = [StandardData; axy(startsWith(string(axy.Identifier), StandardID), :)];
    RefgasData = [RefgasData; axr];  % Reference gas is already extracted by findQ
end

%% --- Store Results in Struct ---
SX.StandardData = StandardData;
SX.RefgasData = RefgasData;
SX.FileNames = FileNames;

%% --- Generate Diagnostic Plots (If Enabled) ---
if strcmp(Plot, 'yes')
    set(0, 'DefaultFigureVisible', 'off');
    fold = "./fig";
    if ~exist(fold, 'dir')
        mkdir(fold)
    end
    subfold = sprintf('./fig/%s',Mode);
    if ~exist(subfold, 'dir')
        mkdir(subfold)
    end

    f = figure('Visible', 'off','Position',[210 55 970 777]);
    residuals = StandardData.RawDeltaValue - StandardData.KnownValue;
    scatter(StandardData.PeakAmplitude, residuals, 60, StandardData.Component, 'filled');
    colorbar;
    xlabel('Amplitude');
    ylabel('Residuals (Measured - Known)')
    saveas(f, sprintf('%s/StandardAmplitudes.png',subfold));

    f = figure('Visible', 'off','Position',[210 55 970 777]);
    tiledlayout(2,2,"TileSpacing",'compact')
    % **Plot 1: Time Series of Standard Measurements**
    nexttile
    hold on;
    uniqueCompounds = unique(StandardData.Component);
    colors = lines(length(uniqueCompounds));
    for i = 1:length(uniqueCompounds)
        comp = uniqueCompounds(i);
        compData = StandardData(StandardData.Component == comp, :);
        plot(compData.Analysis, compData.RawDeltaValue, '-o', 'Color', colors(i,:), 'DisplayName', sprintf('C%d', comp));
    end
    hold off;
    xlabel('Analysis Number');
    ylabel('Raw δ Values (‰)');
    legend('Location', 'best');

    % **Plot 2: Residuals vs. Known Values**
    nexttile
    residuals = StandardData.RawDeltaValue - StandardData.KnownValue;
    scatter(StandardData.KnownValue, residuals, 60, StandardData.Component, 'filled');
    colorbar;
    xlabel('Known δ Values (‰)');
    ylabel('Residuals (Measured - Known)');

    % **Plot 3: Histogram of Standard Deviations**
    nexttile
    histogram(residuals, 'FaceColor', 'blue');
    xlabel('Residual δ (‰)');
    ylabel('Frequency');

    % **Plot 4: Reference Gas Drift Over Time**
    nexttile
    plot(RefgasData.Analysis, RefgasData.RawDeltaValue, '-o', 'Color', 'red', 'LineWidth', 1.5);
    xlabel('Analysis Number');
    ylabel('Reference Gas δ (‰)');

    saveas(f, sprintf('%s/DiagnosticFigs.png',subfold));
    set(0, 'DefaultFigureVisible', 'on');

elseif strcmp(Plot,'no')
end

%% --- Finalize Output ---
SX.Function = 'standardsQ';
end