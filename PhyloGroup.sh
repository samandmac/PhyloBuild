#!/bin/bash

#### START OF THE SCRIPT ####

#It's set up to take some arguements, -s (or --genes) for genes of interest path, -a (or --genomes) for genomes of interest path, -m (or --output) for the output path, and -u (or --clades) for indicating whether to keep or remove clades in the pipeline. You can run without parameters and by default it will use set options.

for arg in "$@"; do #This just sets it so that you can use the double --xxxx and it will convert those to -x.
  shift
  case "$arg" in
    "--genes_interest") set -- "$@" "-p" ;;
    "--genomes_interest") set -- "$@" "-h" ;;
    "--output")   set -- "$@" "-y" ;;
    "--tree_genes")   set -- "$@" "-l" ;;
    "--tree_genomes")   set -- "$@" "-o" ;;
    "--genomes_oi")   set -- "$@" "-t" ;;
    "--genes_oi")   set -- "$@" "-r" ;;
    *)        set -- "$@" "$arg"
  esac
done

#This simply takes the arguement from the user and stores that in a variable.
while getopts "p:h:y:l:o:t:r:" opt
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
		l)
			tree_genes="${OPTARG}" 
			;;
		o)
			tree_genomes="${OPTARG}" 
			;;
		t)
			genomes_oi="${OPTARG}" 
			;;
		r)
			genes_oi="${OPTARG}" 
			;;

	esac
done

shift $((OPTIND-1))

working_directory=`pwd`
genomeInterest=${genome_dir:-${working_directory}/Genomes_Interest}
geneInterest=${gene_dir:-${working_directory}/Carriage_Genes}
genesForTree=${tree_genes:-${working_directory}/Tree_Genes}
genomesForTree=${tree_genomes:-${working_directory}/Tree_Genomes}
plots=${output_dir:-${working_directory}/tmp_output}
genomes_oi=${genomes_oi:-"no"}
genes_oi=${genes_oi:-"no"}

#This checks whether input files are in the correct format or not
python3 py/CheckFileFormat.py $genomesForTree $genesForTree $geneInterest $genomeInterest $genes_oi $genomes_oi 

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

if [ $genomes_oi == "yes" ]
then
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
fi

#Change back to original directory
cd $working_directory

if [ $genomes_oi == "yes" ]
then
	#Copying the genomes that were added in to the same genome folder for making tree.
	echo "Moving genomes added in to template genomes folder..."
	cp -a $genomeInterest/*.fasta $genomesForTree

#We make a list of the genomes we have added in -this is used to both remove the genomes from the default genome folder (in case new genomes are added, the old ones that were tested are removed) and for R plot.
	cd $genomeInterest 
	echo "Making list of genomes added in..."
	ls *.fasta | sed 's/.fasta//g' > $plots/genomes.added.txt
fi

#Change to place where genomes are stored and make a list of genome file names
cd $genomesForTree
echo "Making list of genomes used to make the tree..."
ls *.fasta | sed 's/.fasta//g' > $plots/listGenomes.txt

#We make a selection of genes that are used to help separate the phylogroups of the genomes of interest, and we do this by a specific method in makeGenes.sh. 

#Essentially, we get the sequences for the genes in like the most popular reference E.coli genome. We then run a blastn with these sequences as a query, against the first genome in the list of our genomes of interest. We only keep those sequences for genes that had a match with this sequence, and repeat against the next genome in the list, again keeping only those that have a match. This method repeats until the final genome, after which we transfer the list of gene sequences that had matches for all the tested genomes, and use this to make the tree - it is supposed to help separate out the genomes by their phylogroup.
#TL:DR - makeGenes.sh gets the genes (from a reference genome) and finds orthologs in all genomes being queried. NOTE: We only use 10 results at the moment to save time.
cd $working_directory

cd ClermonTyping/

#We now run a ClermonTyping without Mash - which would also assign phylogroups, but isn't necessary because we don't use it.
for x in `cat $plots/listGenomes.txt`
do
	
	makeblastdb -in $genomesForTree/$x.fasta -input_type fasta -out my_fasta -dbtype nucl #This makes blastdb iteratively for each genome
	
	blastn -query ./data/primers.fasta -perc_identity 90 -task blastn -outfmt 5 -db my_fasta -out my_fasta.xml #This runs blastn querying the primers for ClermonTyping
	
	bin/clermont.py -x my_fasta.xml -s 2000 > temp.txt #This makes a temp.txt file which indicates phylogroup assignment
	
	#Make a variable of the output where we adjust it to the column containing phylogroup assignment
	phylogroup=`head temp.txt | awk '{print $(NF)}'`
	
	#Because albertii, fergusonii are too long for the plot, I've shortened the names.
	if [[ $phylogroup = "albertii" ]]
	then
		phylogroup="ALB"
	elif [[ $phylogroup == "fergusonii" ]]
	then 
		phylogroup="FER"
	elif [[ $phylogroup == "cladeI" ]]
	then 
		phylogroup="I"
	elif [[ $phylogroup == "cladeII" ]]
	then 
		phylogroup="II"
	elif [[ $phylogroup == "cladeIII" ]]
	then 
		phylogroup="III"
	elif [[ $phylogroup == "cladeIV" ]]
	then 
		phylogroup="IV"
	elif [[ $phylogroup == "cladeV" ]]
	then 
		phylogroup="V"
	else
		:
	fi
	#And add the variable to a file, adding the genome name on the same line and storing in plots_extra
	echo -e "${phylogroup}\t$x" >> $working_directory/phylogeny_list.txt
done

cd $working_directory

#Finally, get the user input to determine whether to remove phylogroups with clade assignment or not

echo "Removing strains assigned clade phylogroup"

#Generate new list of genomes to use
python3 py/removeCladesFromList.py $working_directory

#Copy new list over the original list
cp $working_directory/new_genome_list.txt $plots/listGenomes.txt
cp $working_directory/new_phylogeny_list.txt $plots/group_list.txt
cp $plots/group_list.txt $working_directory/group_list.txt

if [ $genomes_oi == "yes" ]
then
	#Removing the genomes of interest from the default genome list in Genomes
	echo "Removing genomes of interest from general genome file"
	cd $genomesForTree
	for f in `cat $plots/genomes.added.txt`
	do 
    		rm ${f}.fasta
	done
fi
cd $working_directory

mv $plots/listGenomes.txt $working_directory/listGenomes.txt
#Removing temporary file
rm -r $plots

rm phylogeny_list.txt
rm new_genome_list.txt
rm new_phylogeny_list.txt
