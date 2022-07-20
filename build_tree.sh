#!/bin/bash

#This is the script where we want to examine the carriage of specific genes of interest over a number of E.coli (strains) genomes.
#To do so, we generate a tree and use this to plot the carriage. The steps of this script are below:

	#1. We need a fasta file of the genes of interest (we've taken this from 1 strain as a reference)
		#For this, I initially used the Clermont scheme genes (i.e. looked at what genes the primers are supposed to select for, and took those). I switched that to a method similar to 		the ISME paper, which iteratively looks for genes present in the genomes and only uses those that are kept in the sequences. The intial genes are all those present in a reference strain.
   	#Therefore, the script initially determines the shared genes in all genome sequences being queried (and those that make up the reference tree!). For the moment, this is only set to the top 10 results (no metric, just the first 10 in the list), but this can be changed as necessary.
		
	#2. We have to determine the phylogrouping of each genome that is being queried (and the references), so therefore the script also runs ClermonTyping and the results of which are saved (see output section).
	
	#3 We run a blastn on each genome (with genes gathered by the ISME method used as queries) to determine which genes are present in the genome and their sequence in the respective genome. The results of this are stored in the BlastResults folder for each genome.
	
	#4. We also require a list of genes of interest as an input file (concatenate the genes of interest into one single txt file, following a specfific formatting i.e:
		#>geneName
		#ATCTAATTATATACATACATATAT etc.
		#Preferably, stick to uncomplicated gene names (no symbols).
	#On this, we will run a blastn to find genomes that have a match to the genes of interest in the list to identify carriage in the tree plot, and then save the list containing both the name of the genes and which genomes have matched with them. The result is a list of genes (per genome) that have been matched to the specific genome. This is saved as a specific output (see output section).
	
	#5. We wish to perform alignment on the blastn results. I use part of Thomas' script for this, which is written in PERL, so I've tried to explain what goes on when we run it. Essentially we put all the results from the blastn (from step 3, step 4 blastn is finished and does not used get further except for the R plot) into one file per gene type, so if a gene is present in the genome file, it gets added to one file of that specific gene type, and then in the next genome if that gene is present the sequence is added to the same file - and so on for each gene in each genome. Essentially, you get a file for each gene and it contains the gene sequences from each genome that had a match. These are saved in the working directory, while they then undergo MAFFT alignment, which generates the files in the FINALTEST.Alignments folder (or whatever you decide to name the directory in the arguments for the perl script) which are, as the name suggests, the aligned gene sequences (of each genome match) for each gene. After this step, all the alignments for all the genes are added to the same file: Final.FINALTEST.aln. This generates a concatenated file of the alignments (note that alignment is before concatenation [of all the genes, not of individual genes though!] which is not exactly what Thomas recommended, but it still works!)
	
	#6. This generates a concatenated file of the alignments, which is great but we aren't done yet. We need to make the tree with PhyML, however we cannot use a .aln file to generate this. So, we need to convert it to a phylip file format, using ALTER. This simply generates a .phy file, which I've designated phylipfor.phy (inventive, I know).
	
	#7. Once we have a .phy file, we can generate a tree using phyml. This generates a whole bunch of files, but the important one is saved to the output folder (it's the one needed for R plots).
	
	#8. Once that's done, we run the R script and generate a plot. Hooray! 

#### START OF THE SCRIPT ####
#Here we begin timing to see how long processes are taking.
start=$SECONDS

#It's set up to take some arguements, -s (or --genes) for genes of interest path, -a (or --genomes) for genomes of interest path, -m (or --output) for the output path, and -u (or --clades) for indicating whether to keep or remove clades in the pipeline. You can run without parameters and by default it will use working_directory/GeneInterest, working_directory/GenomeInterest, working_directory/output, and "yes", respectively. 
for arg in "$@"; do #This just sets it so that you can use the double --xxxx and it will convert those to -x.
  shift
  case "$arg" in
    "--genes") set -- "$@" "-s" ;;
    "--genomes") set -- "$@" "-a" ;;
    "--output")   set -- "$@" "-m" ;;
    "--clades")   set -- "$@" "-u" ;;
    "--labels")   set -- "$@" "-e" ;;
    "--rename_genomes")   set -- "$@" "-l" ;;
    *)        set -- "$@" "$arg"
  esac
done

#This simply takes the arguement from the user and stores that in a variable.
while getopts "s:a:m:u:e:l:" opt
do
	case "$opt" in
		s)
			gene_dir="${OPTARG}"
			;;
		a)
			genome_dir="${OPTARG}"
			;;
		m)
			output_dir="${OPTARG}"
			;;
		u)	
			clades_or_not="${OPTARG}" #added this in, so that indication of keeping clades or not is up to user choice.
			;;
		e)
			labels_or_not="${OPTARG}" #added this in too! It lets user decide if labels should be on the plot or not.
			;;
		l)
			rename_genomes="${OPTARG}" #added this in too! Allows user to rename genomes they added in.
			;;
	esac
done

shift $((OPTIND-1))

#Initialising a few paths as variables
#removed Documents directory variable, as script should work from working directory.
#Made Genomes directory containing genomes for building the tree - these are not user input ones but a template list.
working_directory=`pwd`
genomeInterest=${genome_dir:-${working_directory}/Genomes}
geneInterest=${gene_dir:-${working_directory}/Genes}
genesForTree=$geneInterest/Template_Genes
genomesForTree=$genomeInterest/Template_Genomes
plots=${output_dir:-${working_directory}/output}
blastResults=$working_directory/BlastResults
option=${clades_or_not:-"no"}
labels=${labels_or_not:-"no"}
rename=${rename_genomes:-"no"}

#Making some directories - if they exist already something may have gone wrong with clearing them last time - script may have been ended earlier. 
mkdir -p $plots
mkdir -p $blastResults

#I also clear the phylogroup txt file everytime the script is run, so that it doesn't get appended each time - just in case it hasn't already been removed
echo "Removing phylogroup.txt file in $plots - if it exists"
rm $plots/phylogeny_list.txt

#And gene match list - again this would be different with different genes/genomes.
echo "Removing gene match list from $plots - if it exists"
rm $plots/geneMatches.txt

#First off, I'm removing parts of the names in the genomes we've been given - they aren't very good names but these are the only parts we can delete and still have unique ID's.
#Changing to the directory containing genomes to be added in / checked.
cd $genomeInterest

#For example, changes EC0_06_19646_7#6.contigs_velvet.fa to EC0_06_19646_7.fasta. We only change extension to .fasta because the genome folder uses these extensions, which makes it easier in the script for the loops.
#Side note, this isn't going to affect anything outside of the particular format that some genomes were given in. If the format changes to one more suitable, that's great.
for filename in *.fa; do 
   mv ${filename} "${filename%#*}"
   mv ${filename%#*} "${filename%#*}.fasta"
   echo "Changed $filename to NAME.fasta"
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

#We also want a txt list of the genes of interest, after making the list of genes we also generate a list of those (see later)
cd $working_directory
echo "Making list of genes of interest..."
sed -n -e '/>/p' $geneInterest/geneOI.txt | sed 's/^.*>//' > $plots/genesOI.txt

#Change to place where genomes are stored and make a list of genome file names
cd $genomesForTree
echo "Making list of genomes used to make the tree..."
ls *.fasta | sed 's/.fasta//g' > $plots/List.genomes.txt

#We make a selection of genes that are used to help separate the phylogroups of the genomes of interest, and we do this by a specific method in makeGenes.sh. 

#Essentially, we get the sequences for the genes in like the most popular reference E.coli genome. We then run a blastn with these sequences as a query, against the first genome in the list of our genomes of interest. We only keep those sequences for genes that had a match with this sequence, and repeat against the next genome in the list, again keeping only those that have a match. This method repeats until the final genome, after which we transfer the list of gene sequences that had matches for all the tested genomes, and use this to make the tree - it is supposed to help separate out the genomes by their phylogroup.
#TL:DR - makeGenes.sh gets the genes (from a reference genome) and finds orthologs in all genomes being queried. NOTE: We only use 10 results at the moment to save time.
echo "=============== Step 1: Making List of Genes ==============="
#WARNING, TEMPORARILY REMOVED ISME METHOD BELOW BECAUSE IT DOESN'T SEEM TO BE GREAT SEPARATION, USING geneList.fasta for now which is the CLermont scheme genes.
#bash $working_directory/blast/makeGenes1.sh -p $plots

#WARNING, below is also necessary for the step 1.
#And here's where we make a list of genes used.
#sed -n -e '/>/p' $genesForTree/geneList.txt | sed 's/^.*>//' > $plots/geneNamesTree.txt

cd $working_directory


#Below, we are running steps 2, 3 & 4. We use the ClermonTyping tool to do this, so in a loop for each genome we change to that directory and run the tool, specfiying the genome file and holding a threshold of 2000 - this is the online tool default, and without this threshold a number of genomes cannot be identified. This generates a genome_phylogroups.txt file, of which we take the phylogroup assignment column and add that to a file with the genome name on the same line. We then move back to the documents directory, run a blastn where the geneList.txt file (which is the genes that were made in makeGenes.sh) are queries, the genomes are the subject, and use coverage of 80% and identity of 70%. For each genome, this generates a fasta format file where the gene name is given, then the matching sequence in the genome. Finally, we run another blastn - this time using the genes of interest (so those that are added in by the user) as a query, and still the genomes are a subject. In this case, we just grab the gene names of those that have a match in the genome. We then take the lines in the file, add the genome name in the same line as the gene name, then put that file in the plot_extras directory.
echo "=============== Steps 2, 3 & 4: Phylogrouping, blastn for tree, blastn for genes of interest ==============="
cd ClermonTyping-master/

#We now run a ClermonTyping without Mash - which would also assign phylogroups, but isn't necessary because we don't use it.
for x in `cat $plots/List.genomes.txt`
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
	echo -e "${phylogroup}\t$x" >> $plots/phylogeny_list.txt
done

cd $working_directory

#Finally, get the user input to determine whether to remove phylogroups with clade assignment or not
if [[ $option = "no" ]]
then
	echo "Removing strains assigned clade phylogroup"
	
	#Delete clades from phylogeny list
	#sed -i '/^I\b/d' $plots/phylogeny_list.txt 
	#sed -i '/^II\b/d' $plots/phylogeny_list.txt 
	#sed -i '/^III\b/d' $plots/phylogeny_list.txt 
	#sed -i '/^IV\b/d' $plots/phylogeny_list.txt 
	#sed -i '/^V\b/d' $plots/phylogeny_list.txt 
	#Generate new list of genomes to use
	#awk '{print $2}' $plots/phylogeny_list.txt | grep -f - $plots/List.genomes.txt > $plots/new_genome_list.txt
	python3 py/removeCladesFromList.py $plots
	#Copy new list over the original list
	cp $plots/new_genome_list.txt $plots/List.genomes.txt
	cp $plots/new_phylogeny_list.txt $plots/phylogeny_list.txt
#If not no (so yes) then break and continue.
else
	:
fi

cd $working_directory

#This is the loop where we can run the steps iteratively for each genome.
for x in `cat $plots/List.genomes.txt`
do

	#We run the blastn for the genes to make the tree and save in BlastResults
	echo "Running blastn on $x to identify gene carriage and sequences for template genes (for the tree)"
	blastn -query $genesForTree/geneList.fasta -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid sseq" | sed 's/^\(.\{0\}\)/\1>/' | tr '\t' '\n' >  $blastResults/$x.fasta 
	
	#This is the blastn just to identify which genomes have a match in the genes that are added in by the user.
	blastn -query $geneInterest/geneOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
	
	#We then add the genome name in every line of the file, to indicate genome matches, and save this in plot_extras
	for f in blastMatch.txt
	do 
		cat blastMatch.txt | sed "s/$/\t$x/" >> $plots/geneMatches.txt 
	done
	
done

#Remove blastMatch.txt, no further use for it.
rm blastMatch.txt

cd $working_directory

#This is step 5, where we run alignment on the gene sequences found in the first blastn. A more detailed description of the script can be found in the steps at the start of the script.
#Seems to be the occasional issue with a gene that has a long sequence only found in a few causes the script to fail, bit more testing might need to be done but it works with what I've had so far. May also be dependent on the genome names, symbols seem to mess around with it (for instance, originally I had the # kept in the genome name, after removal I had no issues, with it I had issues).

#Couple of things, the alignment file seems to differ everytime it's made - I suspect that the order of the genes is changed around each time, because the same sequences are there, but it seems to be differenly ordered. This may change the tree slightly each time (indeed, there may be minor differences in trees that are made the same way, but it just seems like the order of branches changes so perhaps not entirely major? Can try to make my own alignment if necessary).
echo "=============== Step 5: Running alignment ==============="
#The alignBlast.pl script was graciously provided by Thomas Otto of the University of Glasgow - it was slightly adjusted for our needs but remains largely the same. For more information on what this script does, see step 5 in the introduction of the script. Perl is a language I am unfamiliar with - if stuff goes wrong with this section I'll be unhelpful.
perl alignBlast.pl BlastResults ALIGNMENTS
find . -name 'ALIGNMENTS.*.fasta' -delete #removes unecessary .fasta files generated from script

#Below is step 6, where we convert the resulting alignment file to a phylip format for tree producing in phyml.
echo "=============== Step 6: Changing .aln format to .phy ==============="
java -jar $working_directory/ALTER/alter-lib/target/ALTER-1.3.4-jar-with-dependencies.jar -cg -i Final.ALIGNMENTS.aln -ia -io Linux -o phylipFor.phy -of PHYLIP -oo Linux -op PhyML

#Step 7, where we run phyml. We use a bootstrap repition of 100.
echo "=============== Step 7: Producing Phylogenetic Tree ==============="
phyml -i phylipFor.phy -b 100

#We copy the resulting newick file format tree data into the plot extras for R.
cp phylipFor.phy_phyml_tree.txt $plots/phylipFor.phy_phyml_tree.txt

#Actually maybe keep this
#Remove the list.genomes.txt, as this is all the genomes used and isn't further necessary.
#rm $plots/List.genomes.txt

echo "=============== Step 8: Generating Plot through R ==============="
Rscript generatePlot.r $plots $option $labels $rename
cp $plots/finalPlot.EMF $plots/finalPlot.$$.EMF
rm $plots/finalPlot.EMF
echo "Rscript ran, output finalPlot.$$.EMF should be in $plots directory"

python3 py/rScriptParameters.py $plots $option $labels $working_directory

#Remove the tmp directory, used in step 5 and no longer needed.
rmdir tmp

cd $genomesForTree
#Removing the genomes of interest from the default genome list in Genomes
echo "Removing genomes of interest from general genome file"
for f in `cat $plots/genomes.added.txt`
do 
    rm ${f}.fasta
done

cd $working_directory

#Clear the blast results directory
echo "Clearing Blast Result directory in case of changes to gene numbers"
rm -r $blastResults

#Clear the alignments file (should stay constant but in case of updates old genes used will still show up without this)
echo "Removing ALIGNMENTS.Alignments file (if it exists)."
rm -r ALIGNMENTS.Alignments
rm Final.ALIGNMENTS.aln

#And the blast folder (let's say you had less genomes than before, it would still try to use the last file in this script!)
echo "Removing uncessary files from blast (if needed)"
cd blast
for i in * 
do
    if ! grep -qxFe "$i" listToKeep.txt #I designate listToKeep.txt as a txt file containing the file names to keep
    then
        echo "Deleting: $i"
        rm "$i"
    fi
done

cd $working_directory
#Deleting extra files generated by programs that are no longer necessary.
echo "Deleting extra files that are unecessary after the script..."
find . -maxdepth 1 -name 'phylipFor*' -delete

end=$SECONDS
runtime=$((end - $start))
mins=$((runtime / 60))
secs=$((runtime % 60))
hour=0
while (( $mins >= 60 ))
do
	hour=$((hour+1))
	mins=$((mins-60))
done

if (( $hour >= 1 ))
then
	echo "Script took ${hour} hour, ${mins} minutes, and ${secs} seconds" 
else
	echo "Script took ${mins} minutes and ${secs} seconds" 
fi
