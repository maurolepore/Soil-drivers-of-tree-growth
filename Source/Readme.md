All the R source code resides in this directory.

There are two key points:
1) The code will only work with the paths correctly set up
2) The code (which generates various input files, rds files etc) will only work with the raw data put in the relevant data directories

To Run all analyses, including generation of rds files, the script "Run analyses.R" should be run.
Run analyses, in turn, calls "All analyses.R"; individual analyses can be run separately from this file.
Sufficient data exist in the data directories of the package (out of the box) to run the scripts:
"Find best model random effects.R"
"Stepwise model selection.R"
"Model averaging.R"
"Legume Growth Analysis.R"
"Correlations Analyses.R"
