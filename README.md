# minimalmarkers
Code for choosing the minimum set of genetic markers needed to differentiate all samples in a genotyping dataset
Usage is ./select_minimal_markers.pl <input_file_name> 

Sample SeqSNP data for Cider apples is included as a test dataset.  To run this, put the perl script and data in the same directory and type:
./select_minimal_markers.pl AppleGenotypes.csv 

The code scales well and has been used to select the minimal SNP markers from a dataset of 35,000 wheat Axiom markers versus over 2000 varieties in a few hours on a MacBook Pro. The Apple example data supplied has 1286 markers and 260 varieties and runs in a few minutes. 

Also supplied is a script for checking the proportion of varieties in the original genotyping file that are resulved by the selected minimal marker set. 
  To run this, give the minimal marker results file as the first command argument and the original genotype file as the second, e.g.
  ./check_results.pl AppleGenotypes_minimal_markers.txt AppleGenotypes.csv
