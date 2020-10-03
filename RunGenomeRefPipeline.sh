source ~/.brew_profile
mkdir 1_RawReads 2_TrimmedReads 3_Alignment
mv *.gz 1_RawReads
cd 1_RawReads
#Quality Trimming
for i in `ls *_R1.fastq.gz|cut -f1 -d "_"`
do
fastp -a AGATCGGAAGAGC -w 16 -x -g -i $i"_R1.fastq.gz" -I $i"_R2.fastq.gz" -o "../2_TrimmedReads/"$i"_Trimmed_R1.fq.gz"  -O "../2_TrimmedReads/"$i"_Trimmed_R2.fq.gz"
done

#Run BWA aligner in folder 3_Alignment

cd 3_Alignment
ln -s ../2_TrimmedReads/*.gz ./
files=`ls *_Trimmed_R1.fq.gz|cut -d "_" -f1`

for i in $files
do

	R1=$i"_Trimmed_R1.fq.gz"
	R2=$i"_Trimmed_R2.fq.gz"
#Alignment
	bwa mem -t 20 /home/data/GATK_bundle_HG38/Homo_sapiens_assembly38.fasta.gz $R1 $R2 > $i".sam"
#Sam to Bam
	sambamba view -f bam -t 14 -p -S -o $i."bam" $i."sam"
#Sort Bam        
	sambamba sort -t 14 $i."bam"

        rm $i".sam"
        rm $i".bam"

#MarkDuplicates
        java -jar /home/data/swdump/picard-tools-1.141/picard.jar MarkDuplicates I=$i".sorted.bam" O=$i".sorted.DupsRmvd.bam" M=$i".metrics" REMOVE_DUPLICATES=true
#Select based on Mapping quality and Pairing        
	sambamba view -t 14 -p -F "mapping_quality >= 20 and proper_pair" -f bam -o $i".final.bam" $i".sorted.DupsRmvd.bam"
#Quality Stats
        /home/data/swdump/qualimap_v2.2.1/qualimap bamqc -bam $i".final.bam" -nt 14 -nw 1000 --java-mem-size=10G

done

#Depth calculation
mosdepth Sacha.final Sacha.final.bam
mosdepth Juliette.final Juliette.final.bam


#RunFreeBayes.sh
#Variant calling using FreeBayes and Variant annotaiton

       ls *.final.bam > Bam.list
        freebayes -f ../Homo_sapiens_assembly38.fasta -L Bam.list > All.vcf
        vcftools --vcf All.vcf --min-meanDP 20 --minQ 20 --out All_filtered --recode

	~/sw/pindel/pindel -f ../../Homo_sapiens_assembly38.fasta -i samples.cfg -c ALL -T 4 -o sv_rmdup
	~/sw/pindel/pindel2vcf -r ../../Homo_sapiens_assembly38.fasta -R GRCh38 -p sv_rmdup_D -d 20160823 -v sv_rmdup_D.vcf

