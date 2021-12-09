#!/usr/bin/perl

$infile = $ARGV[0];
chomp $infile;

open(VCF, "$infile");

open(OUT, ">$infile.genotypes");

while($line = <VCF>)

{

chomp $line;

($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @data) = split(/\t/, $line);

if($chr =~ /^#CHR/){$start =1;}
if($start >0){$i++;}
 
 if($i ==1)
  {
  ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @data) = split(/\t/, $line);
  $head = join("\t", @data);
  print OUT "Marker\t$head\n";
  } 

if($i >1)
{
print OUT "${chr}_$pos";
foreach $cell(@data)
 {
 @fields = split(/:/, $cell);
 if($fields[0] =~ /^0[\/\|]1/){$het++; print OUT "\t1";}
 elsif($fields[0] =~ /^0[\/\|]0/){$hom++; $call1++; print OUT "\t0";}
 elsif($fields[0] =~ /^1[\/\|]1/){$hom++; $call2++;  print OUT "\t2";}
 else{print OUT "\t-1";}
 }
print OUT "\n";
}

}
