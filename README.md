# minimalmarkers
Code for choosing the minimum set of genetic markers needed to differentiate all samples in a genotyping dataset
Usage is ./select_minimal_markers.pl <input_file_name> <optional map file>

Sample SeqSNP data for Cider apples is included as a test dataset.  To run this, put the perl script and data in the same directory and type:
./select_minimal_markers.pl AppleGenotypes.csv Map_locations.txt

The code scales well and has been used to select the minimal SNP markers from a dataset of 35,000 wheat Axiom markers versus over 2000 varieties in a few hours on a MacBook Pro. 

