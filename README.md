# minimalmarkers
Code for choosing the minimum set of genetic markers needed to differentiate all samples in a genotyping dataset
Usage is ./select_minimal_markers.pl <input_file_name> 

Sample SeqSNP data for Cider apples is included as a test dataset.  To run this, put the perl script and data in the same directory and type:

./select_minimal_markers.pl AppleGenotypes.csv 

The script will first pick the marker which discriminates the maximum number of varieties (iteration 1).
For iteration 2, it will find the marker which discriminates the maximum varieties which were NOT discriminated by the marker in iteration 1.
This process continues until either: 1) There are no markers left, 2) adding more markers doesn't add additional varietal discrimination or 3) all varieties are discriminated.  Output is written to tab delimited text files which can be viewed in Excel or similar.

The Apple example data supplied has 1286 markers and 260 varieties and runs in a few minutes but larger datasets may take much longer!
The example data file is comma separated and uses 0 for AA, 1 for AB and 2 for BB.  It will also accept tab separated data and A, AB, B formatted genotype calls as input.


A more detailed explanation of the approach can be found in our paper: https://doi.org/10.1371/journal.pone.0242940



Also supplied is a script for checking the proportion of varieties in the original genotyping file that are resolved by the selected minimal marker set. 
  To run this, give the minimal marker results file as the first command argument and the original genotype file as the second, e.g.
  ./check_results.pl AppleGenotypes_minimal_markers.txt AppleGenotypes.csv
For the example data, you should find that 23 SNPs will discriminate all varieties except "Willy" versus "Connie" - which are not distinguishable.


