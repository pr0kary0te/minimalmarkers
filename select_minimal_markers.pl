#!/usr/bin/perl

$infile = $ARGV[0];
chomp $infile;

if($infile !~ /[\d\D]/){ die "\n\nNo input file specified\nUsage: ./select_minimal_markers.pl marker_data.txt\n\n";}

#The input data should be tab or comma separated. The first column is the marker name and the subsequent ones are genotyping scores for the varieties, 
#The first row is assumed to be a header with your variety names.


$maxmarkers =1000000000000;
#This is normally set to more than the number of markers in the input file (~ 35000) but if it is set lower then markers are prioritised.
#E.g. if set to 5000, then the top 5000 markers by MAF score will be used and the rest ignored. You will get a warning if this limit is reached.
#This feature is intended to speed up the runtime for really big datasets such as 800K Axiom genotyping.

$min_maf = 0.001;
#MAF is minor allele frequency.  This is set to a low level to include as many markers as possible but exclude rare error calls.
#It probably needs optimising for your data.



$min_call_rate = 0.9;
#Ignore markers with less than this proportion of valid (0, 1 of 2 ) calls. (Bad calls might be encoded as -1 or NaN)


#Initiate some hashes.
my %matrix = ();
my %testmatrix =();


#Open the input data file handle - the input file having been specified on the command line
open(IN, "$infile");
$infile =~ s/\..*/_minimal_markers.txt/;
open(OUT, ">$infile");
$infile =~ s/_minimal_markers.txt/_selected_markers_with_headers.txt/;

open(FULL, ">$infile");

print OUT "CumulativeResolved\tMarkerID\tGenotypePattern\n";

#Parse the file header before going through all the data lines
#The file format for the header is "Marker name\tVariety1 name\tVariety2 name\tVariety3 name...."
$head = <IN>;
chomp $head;
#THis REGEX will cope with either tab or comma separated data but bad things may happen if you have a mixture of these characters present!
($id, @header) = split(/[\t\,]/, $head);

#Check to see if the last cell in the input file header has the label "type" - if it does, then remove this from the header and also from the data below
$lastcell = $header[$hlen -1];
if($lastcell =~ /type/){$last = pop(@header);} 

$header = join("\t", @header);
$hlen = @header;

print FULL "ID\t$header\n";

#Start reading the data here
#The file format for the data is "Marker name -> Variety1 score  -> Variety2 score-> Variety3 score...."
while(<IN>)
{
chomp;
#Again, the REGEX will cope with comma or tab separation
($id, @data) = split(/[\t\,]/);

#Again, check if the last column of the header said "type" and get rid of the last data column if it did. 
if($lastcell =~ /type/){$last = pop(@data);}
%alleles = ();


#Common formats are A, B for homozygotes (AA and BB) and AB for heterozygotes or 0,1,2 for hom AA, het AB and hom BB
#We're going to convert A AB B to 0,1,2 format from here on
foreach $cell(@data)
  {
  #Convert A, AB, B format to 0, 1, 2
  $cell =~ s/^AB$/1/;
  $cell =~ s/^A$/0/;
  $cell =~ s/^B$/2/;

  #Replace any cells which aren't 0, 1 or 2 with "x" - This will convert typical bad/missing data values such as "-1" or "NaN" to "x" 
  #This is imporatant as the data are converted to a single string for each marker so there must be exactly one chracter per column.
  if($cell !~ /^[012]$/) {$cell = "x";} 
  
  #Make a hash list of alleles observed in this current row  which aren't bad "x" calls
  if($cell !~ /x/){$alleles{$cell}++;}
  }

$thislen = @data;
#Check that the header and each data row have the same number of cells, in case of file corruption. Die if not.
if($hlen != $thislen){print "$id has length of $thislen which doesn't match header ($hlen) - check your input file\n"; die;}

#Join up the genotype data cells into one single string of concatenated 0,1,2 
$data = join("", @data);

$n_alleles = keys %alleles;

$fails = 0;

#Count the failed calls so we can work out the good call rate
while($data =~ /x/g){$fails++;}
$callrate = ($thislen - $fails)/$thislen;


#Check the call rate is above threshold and we have more than one allele observed (i.e. polymorphic) 
# and add qualifying IDs to the %pattern2id hash. Note that if any SNPs generate identical call patterns 
# then previous IDs will get overwritten as we work through the rows: this is intentional to prevent time 
# wasting with testing SNPs with identical call patterns.


if($callrate > $min_call_rate && $n_alleles >1)
  {
   $pattern2id{$data} = $id;
  }
}

$n = keys %pattern2id;

print "Loaded marker data for $n distinct patterns\n";
close IN;
#All done with reading data from the SNP input file now.


#Now loop over the distinct SNP patterns to organise them by Minor Allele Frequency (MAF) score. This part is only really relevant if
#$maxmarkers is set to fewer than the actual number of input markers, othewise all get used anyway regardless of MAF ordering. 

# print "ID\tmaf\tcall=0\tcall=1\tcall=2\n";
foreach $pattern(keys %pattern2id)
 {
 #Initialise this hash each iteration time to empty it.
 %charcount = ();
 $id = $pattern2id{$pattern};
 while($pattern =~ /([012x])/g) {$charcount{$1}++;}
 $zero = $charcount{0}; 
 $one = $charcount{1};
 $two = $charcount{2};
 $count = $zero+$one+$two;
 #Logic steps to work out which is the second most common call, which we'll define as the minor allele.  
 if($one >= $zero && $zero >= $two){$minor = $zero;} 
 elsif ($zero >= $one && $one >= $two){$minor = $one;}
 elsif ($zero >= $two && $two >= $one){$minor = $two;}
  
 $maf = $minor/$hlen;
#print "$id\t$maf\t$zero\t$one\t$two\n";
  if ($maf > $min_maf)
  {
  push @{$orderbymaf{$maf}}, $pattern; 
  $pattern2maf{$pattern} = $maf;
  }
}
print "Sorting by minor allele frequency..\n";


$count = 0;
foreach $mafscore(sort {$b<=>$a} keys %orderbymaf)
{
$ref = $orderbymaf{$mafscore};
@ary = @$ref;
$l = @ary;
if($count < $maxmarkers)
 {
 foreach $pattern(@ary){$id = $pattern2id{$pattern}; if($count < $maxmarkers){$selected{$pattern} = $id;} 
   else{print "Max marker count of $maxmarkers reached, ignoring further markers\n";}$count++; }
 }
}



$n = keys %selected;
print "$n distinct SNP patterns selected for constructing the optimal  dataset\n";


#set up n x n matrix for data
my $l = @data;
$x = 0; $y = 0;
while($x < $l){$y = 0; while($y < $l){$matrix{$x}{$y} = 0; $y++;} $x++;}
print "Built initial $l x $l scoring matrix\n";


$cumulative = 0;
#Target size it the number of cells in the top right (or bottom left) half of the matrix - i.e. the number of unique varietal comparisons.
$target_size = ($l * ($l-1))/2;
print "\n$target_size varietal comparisons\n";
$currentscore = 1;

print "\nIteration\tCumulativeResolved\tProportion\tMarkerID\tPattern\n";


#This is the main loop where we iterate through all of the available rows of SNP data and find the one that adds the most new "1s" to the overall scoring matrix
#This loop will exit when the $currentscore value is zero - i.e. adding another row doesn't add anything to the overall matrix score 
while($currentscore > 0)
{
#This hash is the current working score matrix - it will be evaluated for this iteration and its contents added to the overal matrix once we have decided which SNP row is best for this iteration
%testmatrix =();
$bestscore = 0;
$iteration++; 

#For this iteration, go through all the markers and see which one differentiates the most varieties not differentiated in previous iterations
#Note that $pattern1 will be the concatenated marker calls for a given marker so will look something like "011122111100002" - for 15 varieties worth of data
foreach $pattern1(sort keys %selected)
   {
   $score = 0;
   $id = $selected{$pattern1};
   %testmatrix =();
   @pattern1 = split(//, $pattern1);
   #print "$id\t$pattern1\n";
   $i = 0; $j = 0;
   $len = @pattern1;
   while($i < $len)
      {
      $ichar = $pattern1[$i];
      $j = $i +1;
      #The logic here is to loop through the genotype string for this marker row and compare each position with every position to the right of it
      # This will fill out one triangular half of the matrix so we end up comparing col1 with col2 then col3 but we don't waste time going back and comparing col2 to col1.
      while($j <$len)
         {
         $jchar = $pattern1[$j];
         #If this cell in the matrix is currently set to sero, i.e. this pair of varieties (i and j) are unresolved, and their genotypes are valid and differnt, then we can set this cell in the test matrix to 1 (= resolved) - otherwise it remains set to zero. 
         if($matrix{$i}{$j} ==0 && $ichar ne "x" && $jchar ne "x" && $ichar != $jchar) 
             {
             #If vars i and j are different, we can set their entry to "1" in the temporary holding matrix and also add 1 to the overall score for the addition of this marker
             $testmatrix{$i}{$j}++; $score++;  
             }  
         $j++;
         }
      $i++;
      }


  #If the current marker is better than others tested, it becomes the new bestscore and its matrix becomes the one to beat!
  if($score >$bestscore)
      {
      $bestscore = $score; 
      %bestmatrix = %testmatrix; 
      $bestpattern = $pattern1;     
      }
   }

      $id = $pattern2id{$bestpattern}; 
      if($bestscore >0)
         {
         $cumulative += $bestscore;
         $proportion_resolved = $cumulative/$target_size;
         $resolved = sprintf("%.6f", $proportion_resolved);
 
         print "$iteration\t$cumulative\t$resolved\t$id\t$bestpattern\n"; 
         print OUT "$cumulative\t$id\t$bestpattern\n"; 
         @tabbed_pattern = split(//,$bestpattern);
         $tabbed_pattern = join("\t", @tabbed_pattern);
         print FULL "$id\t$tabbed_pattern\n";
         }
        $currentscore = $bestscore;
       $i = 0; $j = 0;
         @pattern1 = split(//, $bestpattern);
       $len = @pattern1;
        
        while($i < $len)
         {
         $j = $i +1;
         while($j <$len)
            {
            $matrix{$i}{$j} += $bestmatrix{$i}{$j};  
            $j++;
            }
         $i++;
         }
$bestscore =0;
%bestmatrix =();


}









