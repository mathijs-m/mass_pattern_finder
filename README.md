# mass_pattern_finder
Scripts to analyze UHPLC-HRMS data and find specific patterns in the data and compare several files.

`mass_pattern_finder.py` looks for masses in mass spectrometry files that correspond to a range of chemical formulas, which can be specified through the command line.
`pattern_comparer.py` can take several output files of `mass_pattern_finder.py` and compare the results within a given mass and time accuracy.

## Requirements
The mass spectrometry files should be saved as .mzXML files to be read by the script. The scripts use the pyteomics package to parse through the .mzXML files and the multiprocessing package to parallelize the calculations.

## How to use
Save your mzXML files in a folder and run the script to find the masses  from the command line using `python mass_pattern_finder.py -elements <element_string> -input <path to input file(s)>` to run `mass_pattern_finder` with the default settings (i.e. searching for chlorination with a 5 ppm mass accuracy). Additional help on command line options can be viewed with `python mass_pattern_finder.py --help`.
