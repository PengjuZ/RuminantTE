#!/bin/bash
module load blast/2.5.0-gcc-4.8.5
RepeatMasker=/BIGDATA2/zju_pjzhao_1/00-Software/RepeatMasker/RepeatMasker
RepeatMaskerLib=/BIGDATA2/zju_pjzhao_1/00-Software/RepeatMasker/Libraries/RepeatMasker.lib
##############################################
List="AfricanBuffalo Argali BarbarySheep BlackMuntjac BlueSheep BlueWildebeest BohorReedbuck Bongo Bushbuck Cattle ChineseMuntjac ChineseWaterDeer CommonDuiker CommonEland DefassaWaterbuck ForestMuskDeer Gemsbok Gerenuk Giraffe Goat GrantGazelle GreaterKudu Hartebeest HarveyDuiker HimalayanMuskDeer Ibex Impala IndianMuntjac KirksDikdik Klipspringer LesserKudu MaxwellsDuiker Milu MountainNyala MouseDeer Okapi Oribi Pronghorn PrzewalskiGazelle Reindeer Roedeer RoyalAntelope Sheep Sitatunga Springbok Steenbok Suni ThomsonGazelle TibetanAntelope Topi WaterBuffalo WhiteLippedDeer WhiteTailedDeer Yak"
Dir=/GPUFS/zju_pjzhao_1/02-Ruminant/
for i in $List
do
mkdir $Dir/01-OUT/$i/
cd /GPUFS/zju_pjzhao_1/00-Software/LongRepMarker_v2.0-master/
java -Xmx170G LongRepMarker -r $Dir/00-Fasta/Cattle.fa -T no -k 49 -m 100 -t 10 -R fast -o $Dir/01-OUT/$i/
cp -r /GPUFS/zju_pjzhao_1/00-Software/LongRepMarker_v2.0-master/results/ $Dir/01-OUT/$i/
echo $i OK! >> $Dir/00-Com.OUT.txt
$RepeatMasker -lib $Dir/$i/results/finalRepeats.fa.classified -pa 24 -nolow -gff -dir $Dir/$i/ /BIGDATA2/zju_pjzhao_1/01-Ruminant-TE/00-Ruminant/$i.fa
blastn -db $RepeatMaskerLib -num_threads 24 -out $Dir/$i/results/finalRepeats.fa.classified.blast -max_target_seqs 1 -max_hsps 1 -outfmt 6 -evalue 0.000001 -query $Dir/$i/results/finalRepeats.fa.classified
awk 'ARGIND==1{c[$1]=$2}ARGIND==2{H=$10"#"$11;print $2,int($2+1),$5,$6,$7,$7-$6,$10,$11,c[H]}' $Dir/$i/results/finalRepeats.fa.classified.blast $Dir/$i/$i.fa.out | sed 's/ /\t/g' | grep "NODE" > $Dir/$i/$i.tab.out
cut -f 9 $Dir/$i/$i.tab.out | sort | uniq > $Dir/$i/$i.list.out
awk 'ARGIND==1{A[$9]+=1;B[$9]+=$6;C+=$6}ARGIND==2{print $1"\t"A[$1]"\t"B[$1]"\t"B[$1]/C}' $Dir/$i/$i.tab.out $Dir/$i/$i.list.out | sort -nr -k3 | awk 'ARGIND==1{B[$9][$2]+=$6}ARGIND==2{print $0,0+B[$1][1],0+B[$1][2],0+B[$1][3],0+B[$1][4],0+B[$1][5],0+B[$1][6],0+B[$1][7],0+B[$1][8],0+B[$1][9],0+B[$1][10],0+B[$1][11],0+B[$1][12],0+B[$1][13],0+B[$1][14],0+B[$1][15],0+B[$1][16],0+B[$1][17],0+B[$1][18],0+B[$1][19],0+B[$1][20],0+B[$1][21],0+B[$1][22],0+B[$1][23],0+B[$1][24],0+B[$1][25],0+B[$1][26],0+B[$1][27],0+B[$1][28],0+B[$1][29],0+B[$1][30],0+B[$1][31],0+B[$1][32],0+B[$1][33],0+B[$1][34],0+B[$1][35],0+B[$1][36],0+B[$1][37],0+B[$1][38],0+B[$1][39],0+B[$1][40],0+B[$1][41],0+B[$1][42],0+B[$1][43],0+B[$1][44],0+B[$1][45],0+B[$1][46],0+B[$1][47],0+B[$1][48],0+B[$1][49],0+B[$1][50]}' $Dir/$i/$i.tab.out - | sed 's/ /\t/g' > $Dir/$i/$i.list.divergence.out
done
##############################################
