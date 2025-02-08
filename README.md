## IsoNQ
MATLAB package for Qtegra-generated data file processing.

## Citation
If you use IsoNQ in your research, please cite the following:

pending ...

## Before Using IsoNQ
1. Download and Set Up the Package
Download the IsoNQ package from GitHub and extract the files to a desired location.
(Optional) Move the directory named functions to a preferred location on your computer.

2. Add IsoNQ to Your MATLAB Path
To ensure MATLAB can access the IsoNQ functions, add the directory to your environment path using:

```Matlab
addpath('path_to_functions_directory');
savepath;
```

3. Set the MATLAB Working Directory
This tutorial is designed to be run from the IsoNQ parent directory.
To ensure all file paths work correctly, set your MATLAB working directory to the IsoNQ installation path:

```Matlab
cd('path_to_IsoNQ_directory');
```

You can also set this manually via MATLAB's file explorer.

4. Export Data from Qtegra
IsoNQ's parseQ function processes CSV files exported directly from Qtegra.
To ensure compatibility, follow these steps:

- Open the Query tab in the Qtegra LabBook for the analysis sequence.
- Select all available checkboxes to include all data.
- Click Refresh to update the table.
- In the toolbar, click Export.
- In the export window, select Qtegra EXPORT CSV as the format.
- Click Export and save the file.
- Do not modify the exported file before using it in IsoNQ.

5. View Function Documentation
To view function documentation within MATLAB, use:

```Matlab
help processQ
```

Replace processQ with any IsoNQ function name to see its specific documentation.

## Loading Data Files
This example uses test datasets containing 20 A6 and 3 B4 analyses per sequence.
B4 is used as the reference standard to correct A6 measurements.

Modify the following file paths to match your dataset:

```Matlab
fnames = ["./data/CSIA_A6_B4_test.csv",
          "./data/CSIA_A6_B4_test2.csv"];
```

## Defining Input Parameters
Input parameters are configured using infoQ.
Modify these settings as needed for your dataset:

```Matlab
[IX] = infoQ('SampleRT', [NaN 861 967 1088 1222 1364 1511 1659 1807 1953 2096 2235 2370 2501 2629 2756 2872 2990 3103 3211], ...
             'RefgasRT', [37.829 87.571 137.31 187.06 4008.4 4058.2 4107.9 4157.6], ...
             'Comp', 16:30, ...
             'AddComp', 31:35, ...
             'Report', [-26.15 -31.88 -32.7 -31.99 -33.97 -28.83 -33.77 -33.37 -32.13 -28.46 -32.94 -30.49 -33.2 -29.1 -29.84], ...
             'ReportSD', [0.02 0.02 0.01 0.01 0.02 0.02 0.02 0.03 0.02 0.02 0.01 0.01 0.01 0.01 0.01], ...
             'StandardID', 'B4', ...
             'DetWin', 5, ...
             'RefgasDetWin', 15, ...
             'SampleIDWin', 3, ...
             'Age', './data/sample_ma.xlsx', ...
             'Weight', [25 27 29 31], ...
             'ErrorThreshold', 1, ...
             'OutputFile', 'WeightedMean_output.xlsx', ...
             'PolyN', 2, ...
             'Mode', 'C', ...
             'RemoveFlags', 'yes');
```

To check the input parameters, use:

```Matlab
help infoQ
```

## Running Full Data Processing
To process all data in one step, use:

```Matlab
PX = processQ(fnames, IX, 'Plot', 'yes');
disp(PX.Descriptions);
```

The final processed data is stored in PX.

## Running Individual Steps Manually
To run individual processing steps instead of processQ, execute:

```Matlab
p = 'no'; r = 'yes';  % Define plot and flag removal settings

EX = parseQ(fnames); 
AX = findQ(EX, IX); 
QX = correctQ(AX, IX, 'Plot', p);
CX = combineQ(QX, 'Plot', p);
AQ = assignQ(CX, IX, 'Plot', p);
UCX = errorQ(AQ, AX, IX, 'Plot', p); 
WT = weightQ(AQ, IX, 'UCX', UCX, 'Plot', p, 'RemoveFlags', r);
```

Running these steps separately allows for customization and debugging.

## Troubleshooting
If you encounter errors:
Ensure your MATLAB working directory is set correctly (cd('path_to_IsoNQ_directory')).
Verify that Qtegra CSV files are exported correctly.
Use help <function> to check function requirements.
If figures are not generating, ensure 'Plot' is set to 'yes' in function calls.

## Final Notes
The IsoNQ workflow is designed for n-alkane δ¹³C and δD analysis but may be adapted for similar datasets.
If you modify processing settings, always verify outputs using PX.Descriptions.
For further inquiries, open an issue on the GitHub repository.

## License
This project is licensed under the MIT License.



