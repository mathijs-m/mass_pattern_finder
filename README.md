# mass_pattern_finder
Scripts to analyze UHPLC-HRMS data and find specific patterns in the data and compare several files.

`mass_pattern_finder.py` looks for masses in mass spectrometry files that correspond to a range of chemical formulas, which can be specified through the command line.
`pattern_comparer.py` can take several output files of `mass_pattern_finder.py` and compare the results within a given mass and time accuracy.

## Requirements
The mass spectrometry files should be saved as .mzXML files to be read by the script. The scripts use the pyteomics package to parse through the .mzXML files and the multiprocessing package to parallelize the calculations.

## How to use
Save your mzXML files in a folder and run the script to find the masses  from the command line using `python mass_pattern_finder.py -elements <element_string> -input <path to input file(s)>` to run `mass_pattern_finder` with the default settings (i.e. searching for chlorination with a 5 ppm mass accuracy). Additional help on command line options can be viewed with `python mass_pattern_finder.py --help`. The `<element string>` is a string that specifies the range of chemical formulas for which the masses are reported. This string is formatted as `[minimum occurence]-[element]-[maximum occurence]` and concatenated with `_`. Thus, to find e.g. masses correlating with molecules with between 10 and 15 carbons, 15 to 25 hydrogens, 4 to 5 oxygens and 1 chlorine, the string should be `10-C-15_15-H-25_4-O-5_1-Cl-1`. Custom masses can be defined with `[name]:[mass]`. E.g. the presence of a single phenyl ring can be entered with '1-phenyl:77.0391-1'.

After analyzing the mzXML files with `mass_pattern_finder.py`, the results of the various runs (e.g. sample and blank) can be compared using `pattern_comparer.py`. By running `python pattern_comparer.py -input <path to input file(s)> -reference <path to reference file(s)> -plot_mass_range [min mass]-[max mass] -plot_time_range [min_time]-[max_time]` the file(s) in `<path to input file(s)` are compared to the files in `<path to reference file(s)>`. If the full path of more than one file matches the input file or references file string, all these files are taken as input or reference file. In this case, the script will only report peaks that are within the specified mass and time accuracy in **all** input files and in **one or more** reference files.
