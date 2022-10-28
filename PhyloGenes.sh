#!/bin/bash

#This is the script were we get the genes to make the template of the phylogeny tree. This was given as genesOI.txt - keep the file this way to make PhyloTree work.

for arg in "$@"; do #This just sets it so that you can use the double --xxxx and it will convert those to -x.
  shift
  case "$arg" in
    "--genes_interest") set -- "$@" "-p" ;;
    "--genomes_interest") set -- "$@" "-h" ;;
    "--output")   set -- "$@" "-y" ;;
    "--tree_genes")   set -- "$@" "-r" ;;
    "--tree_genomes")   set -- "$@" "-e" ;;
  esac
done

#See important note below
echo "IMPORTANT!!!!:User should use reference genome called geneList1.txt - this should act as a reference genome which consists of all genes present in the reference genome. For example, we took the nucleotide sequences for genes present in MG1655 for E.coli and align this to the next genome sequence in the genomes file, keeping the genes that exist in both, and aligning this to the next genome sequence, which continues till the final genome is aligned. This can be easily downloaded from NCBI, search for reference strain, get complete sequence option, click \"Send to\" on the overview page for this strain, select \"Gene Features\" in the drop-down menu, choose \"FASTA Nucleotide\" and then press \"Create File\". Once you've obtained this, you can simply transfer the sequence file into your PhyloGenes folder, rename to geneList1.txt, and continue PhyloGene as normal."
sleep 1

echo "===== PhyloGenes will work to determine orthologs in the genomes that can be used to group the genomes for the tree ====="

sleep 2

#This simply takes the arguement from the user and stores that in a variable.
while getopts "p:h:y:r:e:" opt
do
	case "$opt" in
		p)
			gene_dir="${OPTARG}" #User can decide gene directory
			;;
		h)
			genome_dir="${OPTARG}" #User can decide genome directory
			;;
		y)
			output_dir="${OPTARG}" #User can decide output directory
			;;
		r)
			tree_genes="${OPTARG}" 
			;;
		e)
			tree_genomes="${OPTARG}" 
			;;
	esac
done

shift $((OPTIND-1))

#Initialising a few paths as variables, working directory, genes of interest, genome of interest, then tree files, output folder, and a PhyloGenes folder.
working_directory=`pwd`
genomeInterest=${genome_dir:-${working_directory}/Genomes_Interest}
geneInterest=${gene_dir:-${working_directory}/Carriage_Genes}
genesForTree=${tree_genes:-${working_directory}/Tree_Genes}
genomesForTree=${tree_genomes:-${working_directory}/Tree_Genomes}
plots=${output_dir:-${working_directory}/output}
phylogenes=$working_directory/PhyloGenes/

#This checks whether input files are in the correct format or not
python3 py/CheckFileFormat.py $genomesForTree $genesForTree $geneInterest $genomeInterest 

exit_status=$?  # store the exit status for later use

# now lets check the exit status and see if python returned a non-zero exit status
if [ $exit_status -ne 0 ]; then
    echo "An input file  is not in FASTA format, please re-examine your input files and make sure they are in FASTA format"
    exit $exit_status  # exit the bash script with the same status
fi
# continue as usual...
echo "All is good, end of format check"


#Making some directories - if they exist already something may have gone wrong with clearing them last time - script may have been ended earlier. 
mkdir -p $plots

#I also clear the phylogroup txt file everytime the script is run, so that it doesn't get appended each time - just in case it hasn't already been removed
echo "Removing phylogroup.txt file in $plots - if it exists"
rm phylogeny_list.txt

#And gene match list - again this would be different with different genes/genomes.
echo "Removing gene match list from $plots - if it exists"
rm $plots/geneMatches.txt

#First off, I'm removing parts of the names in the genomes we've been given - they aren't very good names but these are the only parts we can delete and still have unique ID's.
#Changing to the directory containing genomes to be added in / checked.
cd $genomeInterest

#For example, changes EC0_06_19646_7#6.contigs_velvet.fa to EC0_06_19646_7.fasta. We only change extension to .fasta because the genome folder uses these extensions, which makes it easier in the script for the loops.
#Side note, this isn't going to affect anything outside of the particular format that some genomes were given in. If the format changes to one more suitable, that's great.
echo "Converting .fa to .fasta - this may not work if the file is already in .fasta - don't worry in this case. In the future, consider using .fasta from the beginning."
for filename in *.fa; do 
   mv "${filename}" "${filename%#*.fa}.fasta"
   #filename=${filename%#*}
   #mv ${filename} "${filename%.fa}.fasta"
   echo "Changed $filename to have .fasta extension"
done

#Change back to original directory
cd $working_directory

#Copying the genomes that were added in to the same genome folder for making tree.
echo "Moving genomes added in to template genomes folder..."
cp -a $genomeInterest/*.fasta $genomesForTree

#We make a list of the genomes we have added in -this is used to both remove the genomes from the default genome folder (in case new genomes are added, the old ones that were tested are removed) and for R plot.
cd $genomeInterest 
echo "Making list of genomes added in..."
ls *.fasta | sed 's/.fasta//g' > $plots/genomes.added.txt

#Change to place where genomes are stored and make a list of genome file names
cd $genomesForTree
echo "Making list of genomes used to make the tree..."
ls *.fasta | sed 's/.fasta//g' > $plots/List.genomes.txt

#Change to that directory, we will work from in here then later transport the orthologs to the Carriage_Genes folder.
cd $phylogenes

echo "Running blast to get some genes which we then use as genes of interest in the main script."

#To convert sequence.txt to a single line for extracting the same sequences from the list of genes to keep
cat geneList1.txt | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' > singleLineGenes.txt

echo "Created singleLinesGenes.txt"

#This sets a number variable which we will use to increment through GeneList$y.txt files - they'll be made as we go!
y=1
numberOfGenes=`ls $genomesForTree | wc -l` 

start=$SECONDS
echo "Beginning to filter genes to those that exist in all genomes"
sleep 1

for x in `cat $plots/List.genomes.txt`
do
	#We're doing this part to say, if it's the last run, then do the below
	if [[ $y == $numberOfGenes ]]
	then
		#This returns percentage identity match, so we can use this to filter our geneList.
		blastn -query geneList$y.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid pident" | sed 's/^\(.\{0\}\)/\1>/' | awk '!seen[$0]++' | sort -k 2n  > Z_Last_File.txt
		#We grab the gene names of the top 10 (lowest percentage match) and get the sequences for these.
		head Z_Last_File.txt | awk '{print $1}' > Z_Ortho_Names.txt
		
		grep -Fwf Z_Ortho_Names.txt singleLineGenes.txt | sed 's/[[:blank:]]*\([^[:blank:]]*\)$/\n\1/' > Z_Orthologs.txt
		
		echo "Extracting the geneList.txt from the orthologs."
	fi
	echo "Running through $x to determine orthologs..."
	#To grab the gene names that appear in the matches -CHANGED TO QCOVS!!!
	blastn -query geneList$y.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid" | sed 's/^\(.\{0\}\)/\1>/' | awk '!seen[$0]++'  > blastGene$y.txt

	
	#To find the same sequences between files AND separate them from one line to multiple lines as in the sequence.txt file
	grep -Fwf blastGene$y.txt singleLineGenes.txt | sed 's/[[:blank:]]*\([^[:blank:]]*\)$/\n\1/' > geneList$(($y + 1)).txt
	
	#Increment y by 1 at the end of the loop
	((y++))

done

#COMMENT
#file=`ls -v | tail -4 | head -1`
#Convert gene list header to normal gene identifier

#COMMENT
sed '/^>/ s/^.*gene=\([Aa-Zz]\+\).*/\1/' Z_Orthologs.txt | sed '1~2s/^/>/' > $genesForTree/geneList.txt

#for only top 10
#COMMENT
#head -20 geneList.txt > $genesForTree/geneList.txt
#mv geneList.txt ~/Documents/Genes/geneList.txt

#COMMENT
#rm $plots/List.genomes.txt

#Removing the genomes of interest from the default genome list in Genomes
#ALL BELOW WAS COMMENTED.
echo "Removing genomes of interest from general genome file"
cd $genomesForTree
for f in `cat $plots/genomes.added.txt`
do 
    rm ${f}.fasta
done

cd $working_directory

rm -r $plots

#And the blast folder (let's say you had less genomes than before, it would still try to use the last file in this script!)
echo "Removing uncessary files from blast (if needed)"
cd $phylogenes
for i in * 
do
    if ! grep -qxFe "$i" listToKeep.txt #I designate listToKeep.txt as a txt file containing the file names to keep
    then
        echo "Deleting: $i"
        rm "$i"
    fi
done

end=$SECONDS
runtime=$((end - $start))
mins=$((runtime / 60))
secs=$((runtime % 60))

echo "It took ${mins} minutes and ${secs} seconds to make the genes for the phylogenetic tree!" 
