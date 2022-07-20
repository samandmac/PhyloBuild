use strict;

##
## for x in `cat List.genomes.txt | head -n 5 ` ; do echo $x; perl in_silico_pcr-master/in_silico_PCR.pl -s ../Genomes/$x.fasta -p primerTab.txt -m 20 -i 10 -r -c 1> output.$x.txt 2> Results2/Res.$x.fasta; done


## for x in `cat ../List.genomes.txt  ` ; do blastn -query ../MSLT.fasta -subject Res.$x.fasta -outfmt 6 -word_size 6 | awk '$9==1' | cut -f 1,2 | perl -nle 'my ($primer,$dir,$seq)=/^(\S+)-(\S+)\s+(\S+)/; if ($dir eq 'P2') { print "$seq\tRevComp" }' > RevComp.$x.txt;  cat Res.$x.fasta | perl ../tree.revComp.pl RevComp.$x.txt > Res.ok.$x.fasta ;done



## read files from directory
## take the genome name $name
## store each amplicon in this specific hash key
## write out the sequences
my %h;
my $dir=shift;
my $res=shift;

my %genomes;

opendir (DIR, $dir) or die "Problem to open opendir $dir: $!\n";

  map {
        if (/(\S+)\.fasta$/){
          my $refName=$1;
	  $genomes{$1}=1;
          print "working on $refName\n";

          # fill the shift hash with the coords
	  open (F, "$dir/$refName.fasta") or die "Problem to read file: $1\n";
	  my $amp;
	  while (<F>){
	      chomp;
	      />(\S+)/;
	      $amp=$1;
	      $_=<F>;
	      $h{$amp}.=">$refName\n$_";
	  }


        }
} readdir(DIR);
closedir(DIR);

system("mkdir $res.Alignments/ ");
system("mkdir tmp");
## Do the alingments and save them in the directory above
foreach my $k (keys %h){
    open (F , "> $res.$k.fasta") or die "Men, you stupid or what: $!\n";
    print F $h{$k};
    close(F);
    !system("mafft --auto $res.$k.fasta > $res.Alignments/$res.$k.aln ") or warn "issues mafft --auto $res.$k.fasta  $!\n";
}


## Combine all the alignments
my %res;


foreach my $k (keys %h){
    ## need to reinitalise the hash every run, 
    my %tmp;
    open (F , "$res.Alignments/$res.$k.aln") or die "Men, I think there is a bug... $res.Alignments/$res.$k.aln  : $!\n";
    my $names;
    my $length;
    while (<F>){
	chomp;
	if (/>(\S+)/){
	    ##  gets the length of the alignement. Assumption that one alignment is set
	    if (defined($names)){ $length=length($tmp{$names}); }
	    $names=$1;
	} else {
	    $tmp{$names}.=$_;
	}
    }
    close(F);

    foreach my $k (keys %genomes){
	if (defined($tmp{$k})){
	    $res{$k}.=$tmp{$k};
	}
	else {
	    $res{$k}.=fillGap($length);
	}
    }
}

## Print final Alignments
open (F , "> Final.$res.aln") or die "Men, you stupid or what: $!\n";
foreach my $k (keys %genomes){
    print F ">$k\n";
    print F $res{$k};
    print F "\n";
}
close(F);

sub fillGap{
    my $l=shift;

    my $res;
    for (1..$l){ $res.="-";}
    return $res;
}

