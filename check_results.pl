#!/usr/bin/perl


open(MARKERSET, "$ARGV[0]");
$headfile = $ARGV[0];

open(FULLDATA, "$ARGV[1]");

while(<MARKERSET>)
{
chomp;
($n, $marker, $data) = split(/[\t\,]/, $_); @data = split(//, $data); $train = @data;
$markers{$marker}++;
}


$head = <FULLDATA>;
chomp $head;
($id, @head) = split(/[\t\,]/, $head);
$full = @head;

while($line = <FULLDATA>)
{
chomp $line;
(@data) = split(/[\t\,]/, $line);
$id = $data[0];
if($markers{$id} >0)
   {
$i++; 
$x = 0; 
$y = 0;
foreach $var1(@head)
      {
      $y = 0;
      $x++;
      foreach $var2(@head)
         {
         $y++;
         $datax = $data[$x];
         $datay = $data[$y];
         if($y > $x){$t++;}
         if($y> $x && $datax =~ /^[012]/ && $datay =~ /^[012]/ && $datax != $datay ){$differs{$var1}{$var2}++;}
         }
      }
   }
}
print "$t comparisons\n";

$x = 0; $y = 0;
foreach $var1(@head)
      {
      $x++;
      $y = 0;
      foreach $var2(@head)
         {
         $y++;
         $paragon = 0;
         if($var1 =~ /paragon/i && $var2 =~ /paragon/i){$paragon =1;}
         if($y> $x && $differs{$var1}{$var2} <1 && $paragon <1 ){$bad{$var1}++; $bad{$var2}++; print "$var1 not resolved from $var2\n";} 
         }
      }

$total = @head;
$n = keys %bad; print "$ARGV[0]\t$n of $total varieties unresolved\n";
