##########################################################################################################
############################## 1. PRE-TREATMENT : FROM FASTQ TO BAM ##################################
#######################################################################################################

#########################################################################################
######################### 1.1. Packages needed to be installed ##########################
#########################################################################################

#For file treatment (bam,bed,bigwig) : samtools, picard tools, bedops, bedtools, deeptols, bedmap
#For peak calling : MACS, SICER, Zerone, HMCan
#For vizualisation : IGV
#Other : bowtie, genome files

#########################################################################################
################################# 1.2. Directories, paths ###############################
#########################################################################################
#Main directory
main='/home/depic2/Documents/Derya/Data'
design='/home/depic2/Documents/Derya/Data/design_A845C148.txt'
cd $main

#Tree creation
mkdir dedup_bam
mkdir bed
mkdir bigwig #Contains the bigwig files we get on internet
mkdir bigwig_bamCompare #Contains bigwig and wig files obtained with bamCompare
mkdir bigwig_bamCoverage #Contains bigwig and wig files obtained with bamCoverage
mkdir -p peaks/Zerone
mkdir -p peaks/HMCan
mkdir csv
mkdir annotation


#Useful paths
PicardDir="/home/depic2/Documents/bioinfo/tools/picard-tools-1.113"
IndexDir="/home/depic2/Documents/bioinfo/Reference_AlignmentGenome_Index/BowtieIndexhg38"
IndexDirColor="/home/depic2/Documents/bioinfo/Reference_AlignmentGenome_Index/BowtieIndexhg38Color"
IndexDirhg38="/home/depic2/Documents/bioinfo/Reference_AlignmentGenome_Index/BowtieIndexhg38"


#HomeMadeTools
Script_Shell="/home/depic2/Documents/bioinfo/Script_Shell"
Script_R="/home/depic2/Documents/bioinfo/Script_R"
Script_Python="/home/depic2/Documents/bioinfo/Script_Python"
patterns=$(cat ~/Documents/bioinfo/Reference_AlignmentGenome_Index/hg38_files/hg38.chrom.sizes | cut -f1)
path_deeptools="/home/depic2/Documents/bioinfo/tools/deepTools-2.5.1/bin"

#peakCallers
path_zerone="/home/depic2/Documents/bioinfo/tools/Zerone"
path_hmcan="/home/depic2/Documents/bioinfo/tools/HMCan"
path_macs="/home/depic2/Documents/bioinfo/tools/MACS"
path_sicer="/home/depic2/Documents/bioinfo/tools/SICER"
#genomeChange
path_liftover="/home/depic2/Documents/bioinfo/tools/liftover"
path_crossmap="/home/depic2/Documents/bioinfo/tools/crossmap"


#Other variables 
genomesize=3049315783 #This is the mapple size of the considered genome : human hg38 here.
line=$(wc -l < $design) #Number of samples
sizes="10000 100000" #Window sizes considered
#Creating lists of chip names, input names and experience names from the design file.
cut -f1 $design > chipList
cut -f2 $design > inputList
cut -f3 $design > expList

#########################################################################################
####################### 1.3. Alignment, sorting and stats ###############################
#########################################################################################

touch BowtieStats.txt
gzip -d *.gz

for n in *.fastq
do 

nom=${n%%.fastq}
echo $nom 
mkdir $nom
bowtie -t -p 10 -S -k 1 -m 1 $IndexDirhg38/hg38 $nom.fastq > $nom/$nom.hg38.sam


##for Solid data use color index and -c option
#bowtie -c -t -p 6 -S -k 1 -m 1 $IndexDirColor/hg38 $nom.fastq > $nom/$nom.sam

x=$(wc -l $n | awk '{print $1}') #nb de lignes dans la fastq

cd $nom

samtools view -@ 6 -b $nom.hg38.sam > $nom.hg38.bam #convert in a smaller file

sambamba_v0.6.5 view -t 10 -f bam -L $hg38Dir/hg38.chrom.bed $nom.hg38.bam > $nom.hg38.sel.bam  ##restrain to classical chromosomes
samtools view -@ 6 -F 4 -b $nom.hg38.sel.bam > $nom.hg38.aligned.bam #ne sélectionne que les alignés

y=4
expr $x / $y >  num_reads #nb reads = x/4 car 4 lignes par reads dans le fichier

sambamba_v0.6.5 view -t 6 $nom.hg38.aligned.bam | wc -l > unique_aligned #nb reads alignés de manière unique
sambamba_v0.6.5 sort -t 6 $nom.hg38.aligned.bam #tri en fonction de la position des reads, créer un fichier aligned.sorted.bam

#remove PCR duplicates

java -jar $PicardDir/MarkDuplicates.jar INPUT=$nom.hg38.aligned.sorted.bam OUTPUT=../dedup_bam/$nom.hg38.sorted.dedup.bam METRICS_FILE=$nom.dupstats REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
sambamba_v0.6.5 index ../dedup_bam/$nom.hg38.sorted.dedup.bam #créer un index *.bai utilisé par IGV
sambamba_v0.6.5 view -t 10 ../dedup_bam/$nom.hg38.sorted.dedup.bam | wc -l > dedu_aligned #nb reads alignés dédupliqués


echo $nom >> ../BowtieStats.txt
cat num_reads >> ../BowtieStats.txt
cat unique_aligned >> ../BowtieStats.txt
cat dedu_aligned >> ../BowtieStats.txt
echo ""  >> ../BowtieStats.txt

rm $nom.hg38.sam
rm $nom.hg38.bam
rm $nom.hg38.sel.bam
rm $nom.hg38.aligned.bam


cd ..

#rm $nom.fastq
done

cd $main/dedup_bam

for file in *bam
do
name=${file%.sorted.dedup.bam}
echo "Converting "$name" from bam to bed"
bam2bed < $file > $main/bed/$name.bed
done





#########################################################################################
##################################### 1.4. Cleaning #####################################
#########################################################################################

for n in */
do
nom=${n%%/*}
mv $n/$nom.sorted.dedup.bam $main/dedup_bam
mv $n/$nom.sorted.dedup.bam.bai $main/dedup_bam
rm $n/*.sam
done


####### Check Bowtie stats and retrieve samples with IP or INPUT below 10 millions reads ############@
