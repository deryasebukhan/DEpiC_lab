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
main="/Users/celinevallot/Documents/bioinfo/DataSets/ChIP-seq/RawData_A806"
cd $main

#Tree creation
mkdir dedup_bam
mkdir bed
mkdir bigwig #Contains the bigwig files we get on internet
mkdir bigwig_bamCompare #Contains bigwig and wig files obtained with bamCompare
mkdir bigwig_bamCoverage #Contains bigwig and wig files obtained with bamCoverage
mkdir peaks
mkdir csv
mkdir annotation


#Useful paths
PicardDir=/Users/celinevallot/Documents/bioinfo/tools/picard-tools-1.113
bamcmpDir=/Users/celinevallot/Documents/bioinfo/tools/bamcmp-master/build
IndexDirhg38=/Users/celinevallot/Documents/bioinfo/Reference_AlignmentGenome_Index/BowtieIndexhg38
IndexDirmm9=/Users/celinevallot/Documents/bioinfo/Reference_AlignmentGenome_Index/BowtieIndexmm9

Script_Shell=/Users/celinevallot/Documents/bioinfo/Scripts/Script_Shell
Script_R=/Users/celinevallot/Documents/bioinfo/Scripts/Script_R
Script_Python=/Users/celinevallot/Documents/bioinfo/Scripts/Script_Python
patterns=$(cat ~/Documents/bioinfo/Reference_AlignmentGenome_Index/hg38_files/hg38.chrom.sizes | cut -f1)
#hg38Dir="/home/depic2/Documents/bioinfo/Reference_AlignmentGenome_Index/hg38_files"
path_deeptools=/Users/celinevallot/Documents/bioinfo/tools/deepTools-2.5.1/bin
path_zerone=/Users/celinevallot/Documents/bioinfo/tools/Rpackage/zerone
path_hmcan=/Users/celinevallot/Documents/bioinfo/tools/HMCan
path_liftover=/Users/celinevallot/Documents/bioinfo/tools/liftover
path_crossmap=/Users/celinevallot/Documents/bioinfo/tools/crossmap

design=/Users/celinevallot/Documents/bioinfo/DataSets/ChIP-seq/design_ChIPseq.txt

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

bowtie -t -p 10 -S -k 1 -m 1 $IndexDirmm9/mm9 $nom.fastq > $nom/$nom.mm9.sam


##for Solid data use color index and -c option
#bowtie -c -t -p 6 -S -k 1 -m 1 $IndexDirColor/hg38 $nom.fastq > $nom/$nom.sam

x=$(wc -l $n | awk '{print $1}')

cd $nom

samtools view -@ 6 -b $nom.hg38.sam > $nom.hg38.bam
samtools view -@ 6 -F 4 -b $nom.hg38.bam > $nom.hg38.aligned.bam
rm $nom.hg38.sam
rm $nom.hg38.bam

samtools view -@ 6 -b $nom.mm9.sam > $nom.mm9.bam
samtools view -@ 6 -F 4 -b $nom.mm9.bam > $nom.mm9.aligned.bam
rm $nom.mm9.sam
rm $nom.mm9.bam


#sambamba_v0.6.5 view -t 10 -f bam -L $hg38Dir/hg38.chrom.bed $nom.bam > $nom.sel.bam ##restrain to classical chromosomes

y=4
expr $x / $y >  num_reads
sambamba_v0.6.5 view -t 6 $nom.hg38.aligned.bam | wc -l > unique_aligned_hg38
sambamba_v0.6.5 view -t 6 $nom.mm9.aligned.bam | wc -l > unique_aligned_mm9
sambamba_v0.6.5 sort -N -t 6 $nom.hg38.aligned.bam
sambamba_v0.6.5 sort -N -t 6 $nom.mm9.aligned.bam

$bamcmpDir/bamcmp -1 $nom.hg38.aligned.sorted.bam -2 $nom.mm9.aligned.sorted.bam -a $nom.hg38.only.bam -A $nom.hg38.first_better.bam -t 6 -n -s mapq

sambamba_v0.6.5 view -t 6 $nom.hg38.only.bam | wc -l > hg38.only
sambamba_v0.6.5 view -t 6 $nom.hg38.first_better.bam | wc -l > hg38.first_better

sambamba_v0.6.5 sort $nom.hg38.only.bam
java -jar $PicardDir/MarkDuplicates.jar INPUT=$nom.hg38.only.sorted.bam OUTPUT=$nom.hg38.dedup.bam METRICS_FILE=$nom.dupstats REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
sambamba_v0.6.5 index $nom.hg38.dedup.bam
sambamba_v0.6.5 view -t 10 $nom.hg38.dedup.bam | wc -l > dedu_aligned
sambamba_v0.6.5 sort $nom.hg38.dedup.bam
mv $nom.hg38.dedup.sorted.bam $main/dedup_bam

bamCoverage -b $nom.hg38.dedup.bam -p 10 -o $main/bwg/$nom.bw -of "bigwig" --normalizeUsingRPKM --extendReads 200 --binSize 50 --blackListFileName $hg38Dir/hg38.blacklist.bed


echo $nom >> ../BowtieStats.txt
cat num_reads >> ../BowtieStats.txt
cat unique_aligned_hg38 >> ../BowtieStats.txt
cat unique_aligned_mm9 >> ../BowtieStats.txt
cat hg38.only >> ../BowtieStats.txt
cat hg38.first_better >> ../BowtieStats.txt
cat dedu_aligned >> ../BowtieStats.txt
echo ""  >> ../BowtieStats.txt

rm $nom.hg38.only.bam
rm $nom.hg38.only.sorted.bam
rm $nom.hg38.dedup.bam
rm $nom.hg38.aligned.bam
rm $nom.hg38.aligned.sorted.bam
rm $nom.mm9.aligned.bam
rm $nom.mm9.aligned.sorted.bam

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
