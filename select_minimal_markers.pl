#!/usr/bin/perl

$infile = $ARGV[0];
chomp $infile;

if($infile !~ /[\d\D]/){ die "\n\nNo input file specified\nUsage: ./select_minimal_markers.pl marker_data.txt\n\n";}

#The input data should be tab or comma separated. The first column is the marker name and the subsequent ones are genotyping scores for the varieties, 
#The first row is assumed to be a header with your variety names.


#Optionally you can supply a file with genetic or physical map locations for your markers as the second command line argument. 
# If you do supply this, the map locations will simply be appended to the ordered marker list created later on in the markers_ordered_at_first_iteration.txt file

$mapfile = $ARGV[1];
if(-e $mapfile)
  {
  open(MAP,"$mapfile"); 
  $head = <MAP>;
  while(<MAP>)
    {
    ($id,$chr,$position,@other) = split(/[\t\,]/, $_);
    $locus = "Chr_$chr $position";
    $id2locus{$id} = $locus;
    }
  }
 


$maxmarkers =1000000000000;
#This is normally set to more than the number of markers in the input file (~ 35000) but if it is set lower then markers are prioritised.
#E.g. if set to 5000, then the top 5000 markers by MAF score will be used and the rest ignored. You will get a warning if this limit is reached.
#This feature is intended to speed up the runtime for really big datasets such as 800K Axiom genotyping.

$min_maf = 0.001;
#MAF is minor allele frequency.  This is set to a low level to include as many markers as possible but exclude rare error calls.
#It probably needs optimising for your data.



$min_call_rate = 0.1;
#Ignore markers with less than this proportion of good (0, 1 of 2 ) calls.  Probably needs optimising.



#Make the file handle hot to prevent write buffering: gives better real-time progress reporting to nohup.
$|=1;


#Initiate some hashes.
my %matrix = ();
my %testmatrix =();


#Open the input data file handle - the input file having been specified on the command line
open(IN, "$infile");
$infile =~ s/\..*/_minimal_markers.txt/;
open(OUT, ">$infile");
print OUT "CumulativeResolved\tMarkerID\tGenotypePattern\n";
open(ORDERED, ">markers_ordered_at_first_iteration.txt");

#Parse the file header before going through all the data lines
#The file format for the header is "Marker name -> Variety1 name-> Variety2 name-> Variety3 name...."
$head = <IN>;
($id, @header) = split(/[\t\,]/, $head);
$hlen = @header;


#Start reading the data here
#The file format for the data is "Marker name -> Variety1 score  -> Variety2 score-> Variety3 score...."
while(<IN>)
{
chomp;
($id, @data) = split(/[\t\,]/, $_); 
%alleles = ();
foreach $cell(@data)
  {
  #Convert A, AB, B format to 0, 1, 2
  $cell =~ s/^AB$/1/;
  $cell =~ s/^A$/0/;
  $cell =~ s/^B$/2/;

  #Replace any cells which aren't 0, 1 or 2 with "x"
  if($cell !~ /^[012]$/) {$cell = "x";} 
  
  #Make a hash list of alleles observed in this current row  which aren't bad "x" calls
  if($cell !~ /x/){$alleles{$cell}++;}
  }

$thislen = @data;
#Check that the header and each data row have the same number of cells, in case of file corruption. Die if not.
if($hlen != $thislen){print "$id has length of $thislen which doesn't match header ($hlen) - check your input file\n"; die;}


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


#Now loop over the distinct SNP patterns to organise them by MAF score. This part is only really relevant if
#$maxmarkers is set to fewer than the actual number of input markers, othewise all get used anyway regardless of MAF ordering. 

 print "ID\tmaf\tcall=0\tcall=1\tcall=2\n";
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
print "$id\t$maf\t$zero\t$one\t$two\n";
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

while($x < $l){$y = 0; while($y < $l){$matrix{$x}{$y} = 0; 
#print "Matrix $l build $x $y\n"; 
$y++;} $x++;
}

print "Built initial $l x $l scoring matrix\n";
$cumulative = 0;

$target_size = ($l * ($l-1))/2;
print "\n$target_size varietal comparisons\n";
$currentscore = 1;

print "\niteration\tCumulativeResolved\tProportion\tMarkerID\tPattern\n";

while($currentscore > 0)
{
%testmatrix =();
$bestscore = 0;
$iteration++; 

#For this iteration, go through all the markers and see which one differentiates the most varieties not differentiated in previous iterations
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
      while($j <$len)
         {
         $jchar = $pattern1[$j];
         if($matrix{$i}{$j} ==0 && $ichar ne "x" && $jchar ne "x" && $ichar != $jchar) 
             {
             $testmatrix{$i}{$j}++; $score++;  
             }  
         $j++;
         }
      $i++;
      }

if($iteration ==1){push @{$orderbyscore{$score}}, $pattern1;}


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

if($iteration == 1)
  {
  foreach $score (sort {$b<=>$a} keys %orderbyscore)
     {
     $ref = $orderbyscore{$score}; @ary = @$ref; 
     foreach $pattern(@ary){$id = $pattern2id{$pattern}; $map = $id2locus{$id}; 
     print ORDERED "$id\t$map\t$score\t$pattern\n";}
     }
  }

}





sub compare
{
my $score = 0;
my $c = $l;
my $x = 0;
my $y = 0;
while($x < $c)
{
$y = 0;
while($y < $c){if($matrix{$x}{$y} ==0 && $testmatrix{$x}{$y} >0){$score++;}$y++; }
$x++;
}
return $score;
}





