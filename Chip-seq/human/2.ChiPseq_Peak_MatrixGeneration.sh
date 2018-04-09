#######################################################################################################
#################### 2. NORMALIZATION FROM BAM FILES AND CREATING COUNT MATRIXES ######################
#######################################################################################################


#### 2.1 CREATING ANNOTATION FILES ("calque que l'on dépose sur les fichiers bam") ####

###### 2.1.1 Fixed windows ######

#For fixed windows method : we take the basic 10k and 100k files and we get rid of blacklist

cd $hg38Dir
bedtools intersect -v -a hg38.10k.bed -b hg38.blacklist.bed > $main/annotation/10k.bed
bedtools intersect -v -a hg38.100k.bed -b hg38.blacklist.bed > $main/annotation/100k.bed

###### 2.1.2 Peak calling ######

########## 2.1.2.1 Zerone ############

#Calling Zerone algorithm and merging 
cd $main/dedup_bam
for n in $(seq 1 $line)
do 
chip=$(head -$n $main/chipList | tail -1)
input=$(head -$n $main/inputList | tail -1)
exp=$(head -$n $main/expList | tail -1)
if [ -f $chip.sorted.dedup.bam ]
then
echo "Doing Zerone for "${exp}
$path_zerone/zerone -w 1000 -l -c 0.95 -0 $main/dedup_bam/$input.sorted.dedup.bam -1 $main/dedup_bam/$chip.sorted.dedup.bam > $main/peaks/Zerone/$exp.1kbp.peak.bed
bedtools merge -d 10000 -i $main/peaks/Zerone/$exp.1kbp.peak.bed > $main/peaks/Zerone/$exp.1kbp.peak_merge.bed
fi
done

#Removing space issues
cd $main/peaks/Zerone
for file in *.bed
do
mv "$file" "${file//[[:space:]]}"
done

########## 2.1.2.2 HMCan ############

#Calling HMCan algorithm (and output log in HMCan folder)
cd $main/dedup_bam
for n in $(seq 1 $line)
do 
chip=$(head -$n $main/chipList | tail -1)
input=$(head -$n $main/inputList | tail -1)
exp=$(head -$n $main/expList | tail -1)
if [ -f $chip.sorted.dedup.bam ]
then
echo "Doing HMCan for "${exp}
HMCan $main/dedup_bam/$chip.sorted.dedup.bam $main/dedup_bam/$input.sorted.dedup.bam $path_hmcan/configurations/config1.txt $exp > $main/peaks/HMCan/log_$exp
fi
done

#Moving HMCan output files (_peaks.narrowPeak and _regions.bed) from dedup_bam to HMCan folder	
mv *.narrowPeak $main/peaks/HMCan
mv *.bed $main/peaks/HMCan

#Removing space issues
cd $main/peaks/HMCan
for file in *.narrowPeak
do
mv "$file" "${file//[[:space:]]}"
done

#Converting .narrowPeak to .bed format
for file in *.narrowPeak
do
name=${file%_*}
cut -f 1-3 $file > $name.peak.bed
done

#Sorting bed files by chr and starting position
for file in *.peak.bed
do
sort-bed $file > sorted_$file
rm $file
rename "s/sorted_//" *
done

#Merging
for file in *peak.bed
do
name=${file%.peak.bed}
bedtools merge -d 10000 -i $file > $name.peak_merge.bed
done

########## 2.1.2.3 Creating the consensus annotation file for Zerone & HMCan ############

#cd $main/annotation
#Rscript $Script_R/Consensus_PeakSet.R
#bedtools merge -d 10000 -i $main/annotation/Zerone_nomerge.bed > $main/annotation/Zerone_other.bed
#bedtools merge -d 10000 -i $main/annotation/HMCan_nomerge.bed > $main/annotation/HMCan_other.bed


#nouveau consensus avec mergeBed
filenames=`find '/home/depic2/Documents/Lou/Data/peaks/Zerone/' -mindepth 1 -maxdepth 1 -type f -print0 | xargs -0 -I {} echo "{}" | grep "peak_merge.bed" `
cat $filenames | sort -k1,1 -k2,2n | mergeBed -i stdin > locations.bed

#pour vérifier que c'est correct:
for f in $filenames
do
awk '{print $0"\t","peakFile-"NR}' $f > ${f}_id
done

filenames_id=`find '/home/depic2/Documents/Lou/Data/peaks/Zerone/' -mindepth 1 -maxdepth 1 -type f -print0 | xargs -0 -I {} echo "{}" | grep "peak_merge.bed_id" `

cat $filenames_id | sort -k1,1 -k2,2n | mergeBed -i stdin - o collapse -c 4 > locations_id.bed


#### 2.2 NORMALIZATION AND CREATING COUNT MATRIXES ####

annot_files="$main/annotation/10k.bed $main/annotation/100k.bed $main/annotation/Zerone_merge.bed $main/annotation/HMCan_merge.bed"



for annot in $annot_files
do
	python $Script_Python/pysam_normalization.py -b $main/dedup_bam -d $main/design_ChIPseq_bam2.txt -a $annot -r $main/csv/Complete
done

#### 2.3 KEEPING ONLY INFO ON SAMPLES (no more IP / Input count infos) ####

Rscript $Script_R keep.R