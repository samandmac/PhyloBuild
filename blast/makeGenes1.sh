#!/bin/bash

#This is the script were we get the genes to make the template of the phylogeny tree.
while getopts "p:" opt
do
	case "$opt" in
		p)
			output_dir="${OPTARG}"
			;;
	esac
done

shift $((OPTIND-1))

#Let's intialise some variables shall we?
docs=~/Documents
genesForTree=$docs/Genes
genomesForTree=$docs/Ecoli/Genomes
plots=$output_dir
blast=$docs/blast

#Change to place where genomes are stored and make a list of genome file names
cd $genomesForTree
ls *.fasta | sed 's/.fasta//g' > $plots/List.genomes.$$.txt

cd $blast
#Running blast to get some genes which we then use as genes of interest in the main script.

#To convert sequence.txt to a single line for extracting the same sequences from the list of genes to keep
cat geneList1.txt | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' > singleLineGenes.txt

#This sets a number variable which we will use to increment through GeneList$y.txt files - they'll be made as we go!
y=1

start=$SECONDS
echo "Beginning to filter genes to those that exist in all genomes!"

for x in `cat $plots/List.genomes.$$.txt`
do

#MISUNDERSTOOD ORIGINAL PAPER, THEY ONLY KEEP THE QUERY GENE THAT HAS 80% COVERAGE AND 70% IDENTITY, SO WE NEED TO APPEND QUERY FILE TO KEEP ONLY THOSE THAT SHOULD BE KEPT
	
	echo "Running through $x to determine orthologs..."
	#To grab the gene names that appear in the matches
	blastn -query geneList$y.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid" | sed 's/^\(.\{0\}\)/\1>/' | awk '!seen[$0]++'  > blastGene$y.txt

	
	#To find the same sequences between files AND separate them from one line to multiple lines as in the sequence.txt file
	grep -Fwf blastGene$y.txt singleLineGenes.txt | sed 's/[[:blank:]]*\([^[:blank:]]*\)$/\n\1/' > geneList$(($y + 1)).txt
	
	#Increment y by 1 at the end of the loop
	((y++))

done

file=`ls -v | tail -4 | head -1`
#Convert gene list header to normal gene identifier
sed '/>/ s/.*_\([A-Z]\{2\}_[0-9]\+.[0-9]\).*/>\1/g' ${file} > geneList.txt
#for only top 10
head -20 geneList.txt > $genesForTree/geneList.txt
#mv geneList.txt ~/Documents/Genes/geneList.txt

rm $plots/List.genomes.$$.txt



end=$SECONDS
runtime=$((end - $start))
mins=$((runtime / 60))
secs=$((runtime % 60))

echo "It took ${mins} minutes and ${secs} seconds to make the genes for the phylogenetic tree!" 

