open(F,"<$ARGV[0]");
while(<F>){chomp;
$i++;
@t=split(/\t/);
$s=$t[3];
$s=$t[4] if $ARGV[0]=~/narrowPeak/;
$s=$t[5] if $ARGV[0]=~/ChromHMM/;
print "$t[0]\t$t[1]\t$t[2]\t$s\t$i\t+\n";
}
close F;
