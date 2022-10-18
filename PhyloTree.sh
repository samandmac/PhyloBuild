#!/bin/bash

#This is the script where we want to examine the carriage of specific genes of interest over a number of E.coli (strains) genomes.
#To do so, we generate a tree and use this to plot the carriage. The steps of this script are below:

	#1. We need a fasta file of the genes of interest (we've taken this from 1 strain as a reference) For this, I initially used the Clermont scheme genes (i.e. looked at what genes the primers are supposed to select for, and took those). I switched that to a method similar to the ISME paper, which iteratively looks for genes present in the genomes and only uses those that are kept in the sequences. The intial genes are all those present in a reference strain. HOWEVER, this has presented issues, and grouping looks pretty great with only the Clermont Scheme - separation of groups is good. Therefore, I am making the ISME Method a separate tool (it isn't quite refined yet) and I will continue with the default being that the user provides a list of genes for tree seperation (by groups, species, or strain) in the form of the geneList.fasta file which contains the nucleotide sequences of the genes used for separation. This should be based off of biological knowledge, or just genes that are good for separation. 
		
	#2. We should also know whether the genomes being queried are separated into specific groups, so therefore we should have a group_list.txt in the working directory which contains the group and then (tab separated) the name of the strain. For instance, using PhyloGroup beforehand for E.coli strains allows the generation of a group_list.txt which contains the phylogroup of strains. This is kept for the plotting script, which uses the group_list.txt file to indicate which strains belong to which group. 
	
	#3 We run a blastn on each genome (with genes gathered by the ISME method used as queries) to determine which genes are present in the genome and their sequence in the respective genome. The results of this are stored in the BlastResults folder for each genome - this will be cleared at the end of the script however.
	
	#4. We also require a list of genes of interest as an input file (concatenate the genes of interest into one single txt file, following a specfific formatting i.e:
		#>geneName
		#ATCTAATTATATACATACATATAT etc.
	#Preferably, stick to uncomplicated gene names (no symbols).
	#On this, we will run a blastn to find genomes that have a match to the genes of interest in the list to identify carriage in the tree plot, and then save the list containing both the name of the genes and which genomes have matched with them. The result is a list of genes (per genome) that have been matched to the specific genome. This is saved as a specific output (see output section).
	
	#5. We wish to perform alignment on the blastn results. I use part of Thomas' script for this, which is written in PERL, so I've tried to explain what goes on when we run it. Essentially we put all the results from the blastn (from step 3, step 4 blastn is finished and does not used get further except for the R plot) into one file per gene type, so if a gene is present in the genome file, it gets added to one file of that specific gene type, and then in the next genome if that gene is present the sequence is added to the same file - and so on for each gene in each genome. Essentially, you get a file for each gene and it contains the gene sequences from each genome that had a match. These are saved in the working directory, while they then undergo MAFFT alignment, which generates the files in the FINALTEST.Alignments folder (or whatever you decide to name the directory in the arguments for the perl script) which are, as the name suggests, the aligned gene sequences (of each genome match) for each gene. After this step, all the alignments for all the genes are added to the same file: Final.FINALTEST.aln. This generates a concatenated file of the alignments (note that alignment is before concatenation [of all the genes, not of individual genes though!] which is not exactly what Thomas recommended, but it still works!)
	
	#6. This generates a concatenated file of the alignments, which is great but we aren't done yet. We need to make the tree with PhyML, however we cannot use a .aln file to generate this. So, we need to convert it to a phylip file format, using ALTER. This simply generates a .phy file, which I've designated phylipfor.phy (inventive, I know).
	
	#7. Once we have a .phy file, we can generate a tree using phyml. This generates a whole bunch of files, but the important one is saved to the output folder (it's the one needed for R plots).
	
	#8. Once that's done, we run the R script and generate a plot. This will be made in the output folder, alongside other data files used in the R script - including the list of gene carriage for the strains, a list of the genes of interest, a list of the genomes used, a list of the genomes that were added to the template folder (and removed), and the Newick format file which contains the tree information.

#### START OF THE SCRIPT ####
#Here we begin timing to see how long processes are taking.
start=$SECONDS

#It's set up to take some arguements, -p (or --genes) for genes of interest path, -h (or --genomes) for genomes of interest path, -y (or --output) for the output path, -l (or --rename_groups) for indicating whether the working directory contains a file called new_names.txt which contains the old names of strains and then new names of strains, for name replacement, -o (or --phylogroup) which is a yes or no setting to indicate whether you used PhyloGroup beforehand, and -t (or --grouping) which indicates whether you included group_list.txt file to provide groups for the strains included in the run. You can run without parameters and by default it will use some defaults - working_directory/Genes, working_directory/Genomes, working_directory/output, "no", "no", and "no", respectively.

for arg in "$@"; do #This just sets it so that you can use the double --xxxx and it will convert those to -x.
  shift
  case "$arg" in
    "--genes_interest") set -- "$@" "-p" ;;
    "--genomes_interest") set -- "$@" "-h" ;;
    "--output")   set -- "$@" "-y" ;;
    "--tree_genes")   set -- "$@" "-l" ;;
    "--tree_genomes")   set -- "$@" "-o" ;;
    "--grouping")   set -- "$@" "-t" ;;
    "--rename_genomes")   set -- "$@" "-r" ;;
    "--phylogroup")   set -- "$@" "-e" ;;
    "--mac")   set -- "$@" "-m" ;;
    "--rename_genes")   set -- "$@" "-g" ;;
    "--phylogenes")   set -- "$@" "-b" ;;
    *)        set -- "$@" "$arg"
  esac
done

#This simply takes the arguement from the user and stores that in a variable.
while getopts "p:h:y:l:o:t:r:e:m:g:b:" opt
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
			grouping_yn="${OPTARG}" 
			;;
		r)
			rename_genomes="${OPTARG}" #added this in too! Allows user to rename genomes they added in.
			;;
		e)
			phylogroup_yn="${OPTARG}"
			;;
		m)
			mac="${OPTARG}" 
			;;
		g)
			rename_genes="${OPTARG}" 
			;;
		b)
			phylogenes="${OPTARG}" 
			;;

	esac
done

shift $((OPTIND-1))

#Initialising a few paths as variables, removed Documents directory variable, as script should work from working directory. Made Genomes directory containing genomes for building the tree - these are not user input ones but a template list.
working_directory=`pwd`
genomeInterest=${genome_dir:-${working_directory}/Genomes_Interest}
geneInterest=${gene_dir:-${working_directory}/Carriage_Genes}
genesForTree=${tree_genes:-${working_directory}/Tree_Genes}
genomesForTree=${tree_genomes:-${working_directory}/Tree_Genomes}
plots=${output_dir:-${working_directory}/output}
blastResults=$working_directory/BlastResults
rename_genomes=${rename_genomes:-"no"}
phylogroup=${phylogroup_yn:-"no"}
grouping=${grouping_yn:-"no"}
mac=${mac:-"no"}
rename_genes=${rename_genes:-"no"}
phylogenes=${phylogenes:-"no"}

#Making some directories - if they exist already something may have gone wrong with clearing them last time - script may have been ended earlier. 
mkdir -p $plots
mkdir -p $blastResults

#I also set it so that if you have used PhyloGroup then you must have a group_list.txt, as this was one of the outputs - so if phylogroup is set to "yes" then we automatically set grouping to "yes" as well.
if [ $phylogroup == "yes" ]
then
	grouping="yes"
	mv $working_directory/List.genomes.txt $plots/List.genomes.txt
fi

#This checks whether input files are in the correct format or not, if files aren't in FASTA format an error will be returned.
python3 py/CheckFileFormat.py $genomesForTree $genesForTree $geneInterest $genomeInterest 

exit_status=$?  # store the exit status for later use

# now lets check the exit status and see if python returned a non-zero exit status, if it does then we'll need to stop
if [ $exit_status -ne 0 ]; then
    echo "An input file  is not in FASTA format, please reexamine your input files and make sure they are in FASTA format"
    exit $exit_status  # exit the bash script with the same status
fi
# continue as usual...
echo "All is good, end of format check"

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
#It's been made a little bit more general - if the gene names are given in .fa format then we switch these to .fasta for convenience.
echo "Converting .fa to .fasta - this may not work if the file is already in .fasta - don't worry in this case. In the future, consider using .fasta from the beginning."
for filename in *.fa; do 
   mv "${filename}" "${filename%#*.fa}"
   mv ${filename} "${filename%.fa}.fasta"
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


#We also want a txt list of the genes of interest, after making the list of genes we also generate a list of those (see later)
cd $working_directory
echo "Making list of genes of interest..."
sed -n -e '/>/p' $geneInterest/genesOI.txt | sed 's/^.*>//' > $plots/genesOI.txt

#Change to place where genomes are stored and make a list of genome file names - if phylogroup was run then this will already exist.
if [ $phylogroup == "no" ]
then
	cd $genomesForTree
	echo "Making list of genomes used to make the tree..."
	ls *.fasta | sed 's/.fasta//g' > $plots/List.genomes.txt
fi

#We make a selection of genes that are used to help separate the phylogroups of the genomes of interest, and we do this by a specific method in makeGenes.sh. 

#Essentially, we get the sequences for the genes in like the most popular reference E.coli genome. We then run a blastn with these sequences as a query, against the first genome in the list of our genomes of interest. We only keep those sequences for genes that had a match with this sequence, and repeat against the next genome in the list, again keeping only those that have a match. This method repeats until the final genome, after which we transfer the list of gene sequences that had matches for all the tested genomes, and use this to make the tree - it is supposed to help separate out the genomes by their phylogroup.
#TL:DR - makeGenes.sh gets the genes (from a reference genome) and finds orthologs in all genomes being queried. NOTE: We only use 10 results at the moment to save time.
#cd $working_directory

#echo "=============== Step 1: Making List of Genes ==============="
#WARNING, TEMPORARILY REMOVED ISME METHOD BELOW BECAUSE IT DOESN'T SEEM TO BE GREAT SEPARATION, USING geneList.fasta for now which is the CLermont scheme genes.
#bash $working_directory/blast/makeGenes1.sh -p $plots

#WARNING, below is also necessary for the step 1.
#And here's where we make a list of genes used.
#sed -n -e '/>/p' $genesForTree/geneList.txt | sed 's/^.*>//' > $plots/geneNamesTree.txt

cd $working_directory


#Below, we are running steps 2 & 3. This is changed to steps 1 & 2 for now. We run a blastn where the geneList.fasta file are queries, the genomes are the subject, and use coverage of 80% and identity of 70%. For each genome, this generates a fasta format file where the gene name is given, then the matching sequence in the genome. Finally, we run another blastn - this time using the genes of interest (so those that are added in by the user) as a query, and still the genomes are a subject. In this case, we just grab the gene names of those that have a match in the genome. We then take the lines in the file, add the genome name in the same line as the gene name, then put that file in the plot_extras directory.
echo "=============== Steps 1 & 2: blastn for tree, blastn for genes of interest ==============="

cd $working_directory

#We need to be aware of which geneList file to use
if [ $phylogenes == "no" ]
then 
	for x in `cat $plots/List.genomes.txt`
	do

			#We run the blastn for the genes to make the tree and save in BlastResults
		echo "Running blastn on $x to identify gene carriage and sequences for template genes (for the tree)"
		blastn -query $genesForTree/geneList.fasta -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid sseq" | sed 's/^\(.\{0\}\)/\1>/' | tr '\t' '\n' >  $blastResults/$x.fasta 
			
		#This is the blastn just to identify which genomes have a match in the genes that are added in by the user.
		blastn -query $geneInterest/genesOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
			
		#We then add the genome name in every line of the file, to indicate genome matches, and save this in plot_extras
		for f in blastMatch.txt
		do 
			cat blastMatch.txt | sed "s/$/\t$x/" >> $plots/geneMatches.txt 
		done
			
	done
fi

#This is the loop where we can run the steps iteratively for each genome.
if [ $phylogenes == "yes" ]
then 
	for x in `cat $plots/List.genomes.txt`
	do

		#We run the blastn for the genes to make the tree and save in BlastResults
		echo "Running blastn on $x to identify gene carriage and sequences for template genes (for the tree)"
		blastn -query $genesForTree/geneList.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid sseq" | sed 's/^\(.\{0\}\)/\1>/' | tr '\t' '\n' >  $blastResults/$x.fasta 
		
		#This is the blastn just to identify which genomes have a match in the genes that are added in by the user.
		blastn -query $geneInterest/genesOI.txt -subject $genomesForTree/$x.fasta -qcov_hsp_perc 80 -perc_identity 70 -outfmt "6 qseqid" > blastMatch.txt
		
		#We then add the genome name in every line of the file, to indicate genome matches, and save this in plot_extras
		for f in blastMatch.txt
		do 
			cat blastMatch.txt | sed "s/$/\t$x/" >> $plots/geneMatches.txt 
		done
		
	done
fi

echo "BLAST has been run over all the genomes."

#Remove blastMatch.txt, no further use for it.
rm blastMatch.txt

cd $working_directory

#This is step 3, where we run alignment on the gene sequences found in the first blastn. A more detailed description of the script can be found in the steps at the start of the script.
#Seems to be the occasional issue with a gene that has a long sequence only found in a few causes the script to fail, bit more testing might need to be done but it works with what I've had so far. May also be dependent on the genome names, symbols seem to mess around with it (for instance, originally I had the # kept in the genome name, after removal I had no issues, with it I had issues). So therefore, it may be that the genes you use for separation need to be somewhat identifiable in all the strains used - however, I use ClermonTyping genes and I have no issues with these - despite some not being present in certain strains. Perhaps this is where biological knowledge should be considered most important.

#Couple of things, the alignment file seems to differ everytime it's made - I suspect that the order of the genes is changed around each time, because the same sequences are there, but it seems to be differenly ordered. This may change the tree slightly each time (indeed, there may be minor differences in trees that are made the same way, but it just seems like the order of branches changes so perhaps not entirely major? Can try to make my own alignment if necessary).
echo "=============== Step 3: Running alignment ==============="
#The align_blast.pl script was graciously provided by Thomas Otto of the University of Glasgow - it was slightly adjusted for our needs but remains largely the same. For more information on what this script does, see step 5 in the introduction of the script. PERL is newer for me, so in the future I'm considering writing a script in Python or bash which will likely have the same functions, albeit hopefully work with any given gene sequences (see above where I complain about the longer gene sequences that don't appear in all strains sometimes causing an issue).
perl align_blast.pl BlastResults ALIGNMENTS
find . -name 'ALIGNMENTS.*.fasta' -delete #removes unecessary .fasta files generated from script

#One slightly annoying thing about this perl script is that concatenation back to the final file for some reason seems to mess with the number of sequences, such that if we ended up using different sequence for the geneList file, then occasionally we get an error which says that the length of sequences differs - which fails the script. Therefore, in addition to the align_blast.pl, I have to run MAFFT on the final concatenated file.
mafft --auto Final.ALIGNMENTS.aln > Final.Aligned.aln

#Below is step 4, where we convert the resulting alignment file to a phylip format for tree producing in phyml.
echo "=============== Step 4: Changing .aln format to .phy ==============="
java -jar $working_directory/ALTER/alter-lib/target/ALTER-1.3.4-jar-with-dependencies.jar -cg -i Final.Aligned.aln -ia -io Linux -o phylipFor.phy -of PHYLIP -oo Linux -op PhyML

#Step 5, where we run phyml. We use a bootstrap repition of 100.
echo "=============== Step 5: Producing Phylogenetic Tree ==============="

if [ $mac == "yes" ]
then
	./PhyML-3.1/PhyML-3.1_macOS-MountainLion -i phylipFor.phy -b 100 --quiet
fi

if [ $mac == "no" ]
then
	phyml -i phylipFor.phy -b 100 --quiet
fi

#We copy the resulting newick file format tree data into the plot extras for R.
cp phylipFor.phy_phyml_tree.txt $plots/phylipFor.phy_phyml_tree.txt

#Step 6 is where the R script is used to generate the plot - in case the user uses the same output folder multiple times I've set it so that the name is changed to a unique number each time. The output will also include the R script with the args switched with the user path - so the user can immediately run the same R script after running the script in case they wish to make any adjustments to the script/plot.
echo "=============== Step 6: Generating Plot through R ==============="
Rscript generate_plot.r $plots $rename_genomes $grouping $rename_genes
cp $plots/finalPlot.EMF $plots/finalPlot.$$.EMF
rm $plots/finalPlot.EMF
echo "Rscript ran, output finalPlot.$$.EMF should be in $plots"

#This is simply the python script which switches args with the user's paths, options, etc.
python3 py/rScriptParameters.py $plots $rename_genomes $working_directory $grouping $rename_genes

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
rm Final.Aligned.aln

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
