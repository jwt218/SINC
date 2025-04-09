%% BEFORE USING ISONQ
% 1. Download the IsoNQ package from GitHub and extract the files to a
%    desired location.
%
% 2. (Optional) Move the directory named 'functions' to a location of your choice on
%    your computer.
%
% 3. In MATLAB, add this directory to the environment path by using:
%    >> addpath('path_to_functions_directory');
%    >> savepath;
%
% 4. This tutorial is designed to be run from the IsoNQ parent directory.
%    To ensure all file paths work as intended, set the current MATLAB
%    working directory to your IsoNQ installation path before running the
%    script. You can do this manually via the MATLAB file explorer or using:
%    >> cd('path_to_IsoNQ_directory');
%
% 5. The parsing tool (parseQ) is designed to process CSV files exported
%    directly from the Qtegra interface without modification. To ensure
%    compatibility:
%
%    - Open the **Query** tab in the Qtegra LabBook for the analysis sequence.
%    - Select **all available checkboxes** to include all data.
%    - Click **Refresh** to update the table.
%    - In the toolbar, click **Export**.
%    - In the export window, select **Qtegra EXPORT CSV** as the format.
%    - Click **Export** and save the file.
%
%    This ensures that parseQ can correctly interpret the file structure without errors.
%
% Use 'help <function>' to see documentation for any function.
% Example: >> help processQ

% The file paths in this tutorial assume that you are running these operations
% from the IsoNQ parent directory. If you are working from a different location,
% update the file paths for 'fnames' and any output directories accordingly.
% Ensure that your MATLAB working directory is set correctly before running.

%% FILE NAMES; QTEGRA EXPORT CSV
% Files contain analysis from sequences with 20 A6 and 3 B4
% each. B4 is used as the reference to correct A6 measurements.
% To apply to your own data files, change these paths.

fnames = ["./data/CSIA_A6_B4_test.csv", ...
    "./data/CSIA_A6_B4_test2.csv"];



%% INPUT PARAMETERS (INFOQ OPTIONS)
% Define key parameters for IsoNQ processing. These include retention times,
% compound selection, reference standards, uncertainty thresholds, and more.
% Modify these parameters as needed for your dataset.

% --- Correction Method Options ---
% 'DriftCorrection' -> 'CompoundSpecific' (default) | 'GlobalDrift'
% 'ScaleCorrection' -> 'Regression' (default) | 'TwoPoint'

% --- Known Values Options ---
% 'KnownValues' -> 'match' (default, follows StandardID)
%                  'A6', 'B4', or 'A7' (built-in standard values)
%                  [custom vector] (user-defined known values)

% --- Outlier Handling Options ---
% 'OutlierMethod' -> 'mad' (default, Median Absolute Deviation)
%                    'zscore' (Z-score outlier removal)
%                    'iqr' (Interquartile Range outlier removal)
%                    'grubbs' (Grubbs' test for extreme values)
%                    'none' (No outlier removal)

% --- Parameter Logging ---
% 'LogParams' -> 'yes' (default) | 'no'
% Saves a timestamped log of all input parameters for reproducibility.

[IX] = infoQ('SampleRT',[NaN 858 964 1085 1219 1361 1508 1656 1804 1950 2093 2231 2367 2498 2625], ...
    'RefgasRT',[37.829 87.571 137.31 187.06 4008.4 4058.2 4107.9 4157.6], ...
    'KnownValues','match', ...
    'StandardID','B4', ...
    'PeakWindow',5, ...
    'RefgasWindow',15, ...
    'SampleMatchLen',3, ...
    'AgeFile','./data/sample_ma.xlsx', ...
    'WeightedComp',[25 27 29], ...
    'UncertaintyThreshold',1.5, ...
    'PolyOrder',2, ...
    'Mode','C', ...
    'RemoveFlags','yes', ...
    'DriftCorrection','CompoundSpecific', ...
    'ScaleCorrection','Regression', ...
    'OutlierMethod','mad', ...
    'OutlierThreshold',5, ...
    'LogParams','yes');

help infoQ



%% RUN ALL PROCESSES
% The final results are stored in PX.

PX = processQ(fnames,IX,'Plot','no');

disp(PX.Descriptions)



%% RUN INDIVIDUAL OPERATIONS (PERFORMED IN PROCESSQ)
% The following functions are executed automatically in processQ.
% Running them individually allows step-by-step inspection of each stage.
% Use 'help <function>' for more details on each operation.

p = 'no'; r = 'yes';

[EX] = parseQ(fnames);
[AX] = findQ(EX,IX);
[QX] = correctQ(AX,IX,'Plot',p);
[CX] = combineQ(QX,'Plot',p);
[AQ] = assignQ(CX,IX,'Plot',p);

[UCX] = errorQ(AQ,AX,IX,'Plot',p);
[WT] = weightQ(AQ,IX,'Error',UCX,'Plot',p,'RemoveFlags',r);



%% WRITE OUTPUT TO FILE
% Write data stored in any IsoNQ output struct to a file of a specified
% directory and format.

writeQ(PX);

writeQ(PX, 'OutputDir', 'Format', 'xlsx', 'Timestamp', false);



%% COMPILE STANDARDS
% Use EX = parseQ(...); for different datasets, then input EX structs to
% standardsQ: e.g., standardsQ(IX,EX1,EX2,EX3).

[EX1] = parseQ(fnames(1));
[EX2] = parseQ(fnames(2));

[SX] = standardsQ(IX,EX1,EX2,'Plot','yes');


%% MERGE MULTIPLE ISONQ DATASETS
% Combine multiple processQ output structs (PX) into one struct. This is
% useful if different standards were analyzed across multiple sessions for
% one profile or group of samples.

PX1 = processQ(fnames(1), IX, 'Plot', 'no');
PX2 = processQ(fnames(2), IX, 'Plot', 'no');
[MX] = mergeQ(PX1, PX2);


%% SORTING DATA WITHIN THE MERGED STRUCT
% To sort tables in the MX struct, try sortrows(...). For example, to sort
% the replicate mean delta values by sample ID, use:

delta_sorted = sortrows(MX.MeanDelta, 'ID');

% Sorting options depend on what column structure is in each table.
% Check available columns with:
disp(MX.MeanDelta.Properties.VariableNames);


%% TROUBLESHOOTING
% If you encounter errors:
% - Ensure you have set the correct working directory.
% - Confirm that Qtegra CSV files are formatted as expected.
% - Use 'help <function>' to check input requirements for any function.
