#!/bin/bash
##############################################
module purge
module load cufflinks/2.2.1-gcc-4.8.5
module load hmmer/3.1b2-icc-14.0.2
module load emboss/6.6.0-gcc-4.8.5
module load blast/2.5.0-gcc-4.8.5
module load diamond/0.9.17
module load bedtools2/2.26.0-gcc-4.8.5
module load samtools/1.11-gcc-4.8.5
Exonerate=Exonerate.pl
Filter=Filter.sh
##############################################
Genome=/BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/00-Ruminant
OUT=/BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/01-ORF
Protein=/BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/01-DB-ZNF/01-CDS.MRD.fa
##############################################
List="AfricanBuffalo Argali BarbarySheep BlackMuntjac BlueSheep BlueWildebeest BohorReedbuck Bongo Bushbuck Cattle ChineseMuntjac ChineseWaterDeer CommonDuiker CommonEland DefassaWaterbuck ForestMuskDeer Gemsbok Gerenuk Giraffe Goat GrantGazelle GreaterKudu Hartebeest HarveyDuiker HimalayanMuskDeer Ibex Impala IndianMuntjac KirksDikdik Klipspringer LesserKudu MaxwellsDuiker Milu MountainNyala MouseDeer Okapi Oribi Pronghorn PrzewalskiGazelle Reindeer Roedeer RoyalAntelope Sheep Sitatunga Springbok Steenbok Suni ThomsonGazelle TibetanAntelope Topi WaterBuffalo WhiteLippedDeer WhiteTailedDeer Yak"

for i in $List
do
samtools faidx $Genome/$i.fa
makeblastdb -dbtype nucl -in $Genome/$i.fa -out $Genome/$i
mkdir $OUT/$i/
blastn -db $Genome/$i -num_threads 24 -out $OUT/$i/$i.blast -outfmt 6 -evalue 0.000001 -query $Protein
cat $OUT/$i/$i.blast | awk '{if($10 < $9){T=$9;$9=$10;$10=T};if($3 > 95){print $1,$2,$9,$10}}' | sort -k1,1 -k2,2 -k3n,3 | awk '{if(($3+$4)/2-SS[$1][$2] < 1000000 && A==$1 && B==$2){SS[$1][$2]=($3+$4)/2;A=$1;B=$2}else{SS[$1][$2]=($3+$4)/2;N+=1;A=$1;B=$2};print $0,N}' > $OUT/$i/$i.S.blast
awk 'ARGIND==1{S[$5]+=0;G[$5]+=0;if(S[$5]=="0"){S[$5]=$3};if(G[$5]<$4){G[$5]=$4}}ARGIND==2{print $0,S[$5],G[$5]}' $OUT/$i/$i.S.blast $OUT/$i/$i.S.blast | awk '{print $1,$2,$6,$7}' | uniq > $OUT/$i/$i.SL.blast
awk 'ARGIND==1{c[$1]=$2}ARGIND==2{if($3-10000 < 0){KS=0}else{KS=$3-10000};if($4+10000 > c[$2]){KE=c[$2]}else{KE=$4+10000};print $0,KS,KE}' $Genome/$i.fa.fai $OUT/$i/$i.SL.blast | sed 's/ /\t/g' > $OUT/$i/$i.SLI.blast
$Exonerate $OUT/$i/$i.SLI.blast $Genome/$i.fa $Protein >> $OUT/$i/$i.SLIR.blast
cat $OUT/$i/$i.SLIR.blast | awk '{if($3 == "gene" || $3 == "exon"){print}}' | sed 's/exonerate:est2genome/Exonerate/g' | sed 's/:/\t/g' | sed 's/-\t/MASKZPJ\t/g' |  sed 's/-/\t/g' | awk -F "\t" '{print $1"\t"$4"\t"$5"\t"$6+$2"\t"$7+$2"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' | sed 's/MASKZPJ\t/-\t/g' | awk '{if($3 == "gene"){ID=$13;print $1,$2,$3,$4,$5,$6,$7,$8,"ID="ID}else{NN[ID]+=1; print $1,$2,$3,$4,$5,$6,$7,$8,"ID=exon"NN[ID]";Parent="ID;print $1,$2,"CDS",$4,$5,$6,$7,$8,"ID=CDS"NN[ID]";Parent="ID}}' | sed 's/ /\t/g' > $OUT/$i/$i.SLIR.gff
gffread $OUT/$i/$i.SLIR.gff -g $Genome/$i.fa -y $OUT/$i/$i.SLIR.Protein.fa
hmmsearch -E 10000000 --domT 0 --incE 10000000 --incdomT 0  --domtblout $OUT/$i/$i.Pfam.txt --cpu 24 /BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/00-Pfam-KZFP/Pfam-KZFP.Yes.hmm $OUT/$i/$i.SLIR.Protein.fa
$Filter $OUT/$i/$i.Pfam.txt | grep -v "#" | awk '{print $1,$4,$11}' | uniq | sed 's/ /\t/g' | sort -k1 > $OUT/$i/$i.Pfam.F.txt
awk 'ARGIND==1{A[$1][$2]=$3}ARGIND==2{A[$1]["KRAB"]+=0;A[$1]["zf-C2H2"]+=0;A[$1]["SCAN"]+=0;A[$1]["DUF3669"]+=0;print $1,A[$1]["KRAB"],A[$1]["zf-C2H2"],A[$1]["SCAN"],A[$1]["DUF3669"]}' $OUT/$i/$i.Pfam.F.txt $OUT/$i/$i.Pfam.F.txt | uniq | awk '{if($2+$4+$5 > 0 && $3 > 0){print}}' | sed 's/ /\t/g' > $OUT/$i/$i.Pfam.Domain.txt
cut -f 1 $OUT/$i/$i.Pfam.Domain.txt | seqtk subseq $OUT/$i/$i.SLIR.Protein.fa - > $OUT/$i/$i.Pfam.Domain.fa
##sixpack -nofirstorf Y -nolastorf Y -sequence $Genome/$i.fa -outseq /BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/01-ORF/$i.ORF6.fa -outfile /BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/01-ORF/$i -orfminsize 50 -mstart Y -nodescription
##getorf -find 1 -minsize 150 -sequence $Genome/$i.fa -outseq /BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/01-ORF/$i.ORF.fa
##hmmscan --domtblout /BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/01-ORF/$i.Pfam.6.txt --cpu 24 /BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/00-Pfam-KZFP/Pfam-KZFP.Yes.hmm /BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/02-KZFP/01-ORF/Cattle.ORF6.fa
echo $i OK!
done

