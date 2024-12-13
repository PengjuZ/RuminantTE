#!/bin/bash
####################################################################################
#conda create -n Kmers
#source /share/apps/anaconda3/bin/activate Kmers
#conda install -c bioconda jellyfish
#conda install -c bioconda kmergenie
source /share/apps/anaconda3/bin/activate Ntools
####################################################################################
Genome=/share/home/zju_zhaopj/02-RuminantTE/00-Genome
TEanno=/share/home/zju_zhaopj/02-RuminantTE/01-TEanno
OUT=/share/home/zju_zhaopj/02-RuminantTE/05-TEinfor

###01-infors
mkdir $OUT/01-infors
while read line
do
awk '{if($16 != "*"){$16="#"};print $5,$6,$7,$9,$10,$11,$15,$16}' $TEanno/$line/RepeatMasker/$line.fa.out | grep "NODE" | sed 's/ /\t/g' > $OUT/01-infors/$line.01.infor.bed
bedtools intersect -a $OUT/01-infors/$line.01.infor.bed -b $OUT/01-infors/$line.01.infor.bed -wo | awk '{if($5 != $13)print}' > $OUT/01-infors/$line.02.compare.bed
awk '{if($8 == "#" && $16 == "*"){a=$1;b=$2;c=$3;d=$4;e=$5;f=$6};if($8 == "*" && $16 == "#"){a=$9;b=$10;c=$11;d=$12;e=$13;f=$14};if($8 == $16 && $3-$2 > $11-$10){a=$9;b=$10;c=$11;d=$12;e=$13;f=$14};if($8 == $16 && $3-$2 < $11-$10){a=$1;b=$2;c=$3;d=$4;e=$5;f=$6};print a,b,c,d,e,f}' $OUT/01-infors/$line.02.compare.bed | sort | uniq > $OUT/01-infors/$line.03.rm.bed
awk 'ARGIND==1{A[$1][$2][$3]=1}ARGIND==2{if(A[$1][$2][$3]!=1);if($3 - $2 >= 100){print}}' $OUT/01-infors/$line.03.rm.bed $OUT/01-infors/$line.01.infor.bed > $OUT/01-infors/$line.03.ok.bed
cut -f 1,2,3,4 $OUT/01-infors/$line.03.ok.bed | sed 's/\tC/\t-/g' | bedtools getfasta -fi $Genome/$line.fa -bed - | seqkit fx2tab | grep -v "NN" | awk '{print ">"$1; print$2}' > $OUT/01-infors/$line.04.RE.fa
RepeatMaskerfa=/share/home/zju_zhaopj/.conda/envs/RepeatMasker/share/RepeatMasker/Libraries/RepeatMasker.lib.f/RepeatMasker.lib
blastn -db $RepeatMaskerfa -query $line.04.RE.fa -outfmt 7 -num_threads 24 -max_target_seqs 1 | grep -v "# " > $line.05.blast.out
cut -f 1,2 $line.05.blast.out | uniq | cut -f 2 | sort | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' | sort -nr -k1 > $line.06.blast.stats.txt
awk 'ARGIND==1{A[$2]+=$4}ARGIND==2{print $2"\t"A[$2]}' $line.05.blast.out $line.06.blast.stats.txt | sort -nr -k2 > $line.06.length.stats.txt
cut -f 1,2 $line.05.blast.out | uniq | cut -f 2 | sed 's/#/\t/g' | cut -f 2 | sort | uniq -c| sed 's/^ *//g' | sed 's/ /\t/g' | sort -nr -k1 > $line.07.family.blast.stats.txt
sed 's/#/\t/g' $line.05.blast.out | awk 'ARGIND==1{A[$3]+=$5}ARGIND==2{print $2"\t"A[$2]}' - $line.07.family.blast.stats.txt | sort -nr -k2 > $line.07.family.length.stats.txt
done < /share/home/zju_zhaopj/02-RuminantTE/03-Tree.list.txt
awk '{if($1 >= 1000)print}' *.07.family.blast.stats.txt | cut -f 2 | sort -r | uniq | grep -v -E "snRNA|Satellite" > $OUT/01-infors.family.txt
#length-stats
while read line
do
Input1=$OUT/01-infors/$line.07.family.length.stats.txt
grep -v -E "snRNA|Satellite" $Input1 > $OUT/01-infors-stats/Family/$line.txt
Input2=$OUT/01-infors/$line.06.length.stats.txt
grep -v -E "snRNA|Satellite" $Input2 > $OUT/01-infors-stats/SubFamily/$line.txt
done < /share/home/zju_zhaopj/02-RuminantTE/03-Tree.list.txt
#
while read line
do
Input1=$OUT/01-infors-stats/Family/$line.txt
awk 'ARGIND==1{S+=$2}ARGIND==2{if($2/S > 0.05){print $1}}' $Input1 $Input1 
done < /share/home/zju_zhaopj/02-RuminantTE/03-Tree.list.txt | sort | uniq > Family.TE.list.txt
while read line
do
Input1=$OUT/01-infors-stats/Family/$line.txt
awk 'ARGIND==1{M1[$1]=1}ARGIND==2{S+=$2;if(M1[$1]==1){M+=$2}}ARGIND==3{if(M1[$1]==1){print $1,$2,S,$2/S;print "Other",S-M,S,(S-M)/S}}' Family.TE.list.txt $Input1 $Input1 | \
     sort | uniq | awk '{print ID,$0}' ID=$line - | sed 's/ /\t/g' | sort -k1,1 -k5nr,5
done < /share/home/zju_zhaopj/02-RuminantTE/03-Tree.list.txt > Family.txt
#
while read line
do
Input1=$OUT/01-infors-stats/SubFamily/$line.txt
awk 'ARGIND==1{S+=$2}ARGIND==2{if($2/S > 0.001){print $1}}' $Input1 $Input1 
done < /share/home/zju_zhaopj/02-RuminantTE/03-Tree.list.txt | sort | uniq > SubFamily.TE.list.txt
while read line
do
Input1=$OUT/01-infors-stats/SubFamily/$line.txt
awk 'ARGIND==1{M1[$1]=1}ARGIND==2{S+=$2}ARGIND==3{if(M1[$1]==1){print $1,$2/S}}' SubFamily.TE.list.txt $Input1 $Input1 | \
     sort | uniq | awk '{print ID,$0}' ID=$line - | sed 's/ /\t/g' | sort -k1,1 -k5nr,5
done < /share/home/zju_zhaopj/02-RuminantTE/03-Tree.list.txt > SubFamily.txt

###02-infors
source /share/apps/anaconda3/bin/activate Ntools
mkdir $OUT/02-infors
while read line1
do
mkdir $OUT/02-infors/$line1
     while read line2
     do
     Family=$(echo $line2 | cut -d / -f 2 )
     grep "$line2" $OUT/01-infors/$line1.05.blast.out | cut -f 1 | sed 's/:/\t/g' | sed 's/-/\t/g'> $OUT/02-infors/$line1/$Family.01.bed
     awk 'ARGIND==1{A[$1][$2][$3]=1}ARGIND==2{if(A[$1][$2][$3] == 1)print}' $OUT/02-infors/$line1/$Family.01.bed $OUT/01-infors/$line1.03.ok.bed > $OUT/02-infors/$line1/$Family.02.ok.bed
     cut -f 1,2,3,4 $OUT/02-infors/$line1/$Family.02.ok.bed | sed 's/\tC/\t-/g' | bedtools getfasta -fi $Genome/$line1.fa -bed - | seqkit fx2tab | grep -v "NN" | awk '{print ">"$1; print$2}' > $OUT/02-infors/$line1/$Family.03.fa
     done < $OUT/01-infors.family.txt
done < /share/home/zju_zhaopj/02-RuminantTE/03-Tree.list.txt

###03-merge
mkdir $OUT/03-infors-merge
while read line1
do
     while read line2
     do
     Family=$(echo $line2 | cut -d / -f 2 )
     Fasta=$OUT/02-infors/$line1/$Family.03.fa
     awk '{H=substr($1,0,1);if(H == ">"){S+=1;print ">"ID"_"S}else{print}}' ID=$line1 $Fasta >> $OUT/03-infors-merge/$Family.TE.fa
     awk '{H=substr($1,0,1);if(H == ">"){S+=1;print $1"\t"ID"_"S}}' ID=$line1 $Fasta >> $OUT/03-infors-merge/$Family.TE.list.txt
     done < $OUT/01-infors.family.txt
done < /share/home/zju_zhaopj/02-RuminantTE/03-Tree.list.txt
###
while read line2
do
     Family=$(echo $line2 | cut -d / -f 2 )
     grep ">" $OUT/03-infors-merge/$Family.TE.fa | sed 's/_/\t/g' | sed 's/>//g' | cut -f 1 | uniq -c | sed 's/^ *//g' | sed 's/ /\t/g' > $OUT/03-infors-merge/$Family.stats.txt
done < $OUT/01-infors.family.txt

###04-infors-kmer
source /share/apps/anaconda3/bin/activate Kmers
mkdir $OUT/04-infors-kmer
while read line2
do
     Family=$(echo $line2 | cut -d / -f 2 )
     mkdir $OUT/04-infors-kmer/$Family
     for var in {5..200}  
     do
     ntcard -t 1 -k $var -p $OUT/04-infors-kmer/$Family/out $OUT/03-infors-merge/$Family.TE.fa
     head -3 $OUT/04-infors-kmer/$Family/out\_k$var.hist | awk '{print ID,$1,$2}' ID=$var - | sed 's/ /\t/g' >> $OUT/04-infors-kmer/$Family.kmer.txt
     done 
done < $OUT/01-infors.family.txt

###05-infors-kmer
mkdir $OUT/05-infors-kmer-run
cd $OUT/05-infors-kmer-run
while read line2
do
     Family=$(echo $line2 | cut -d / -f 2 )
     grep "F0" $OUT/04-infors-kmer/$Family.kmer.txt | sort -nr -k3 > $OUT/04-infors-kmer/$Family.kmer.sort.txt
     K=$(cat $OUT/04-infors-kmer/$Family.kmer.sort.txt | head -1 | cut -f 1)
     echo $Family $K
     jellyfish count -C -m $K -s 1G -c 7 -t 12 -o $Family.jf $OUT/03-infors-merge/$Family.TE.fa 
     jellyfish dump -c -t $Family.jf > $Family.tsv
     jellyfish stats $Family.jf -o $Family.stats.txt
     jellyfish histo -t 12 $Family.jf > $Family.histo
     T=$(cut -f 2 $Family.tsv | sort -nr | head -1000 | tail -1)
     awk '{if($2 >= T)print}' T=$T $Family.tsv | sort -nr -k2 > $Family.top.tsv
     awk '{print ">id"NR;print $1}' $Family.top.tsv > $Family.top.fa
     samtools faidx $Family.top.fa
     rm -rf $Family.tsv
     rm -rf $Family.jf
     bwa index -b 500000000 $OUT/03-infors-merge/$Family.TE.fa
     bwa aln -n 0 -o 0 -R 50000000 -t 8 -f $Family.top.sai $OUT/03-infors-merge/$Family.TE.fa $Family.top.fa
     bwa samse -n 50000000 -f $Family.top.sam $OUT/03-infors-merge/$Family.TE.fa $Family.top.sai $Family.top.fa  
     grep -v "SN:" $Family.top.sam | grep -v "@PG" | sed 's/X0:i://g' | sed 's/XA:Z://g' | awk '{print $1,$10,$14,$3";"$20}' | tr ';' '\n' > $Family.tab.sam
     awk '{if($2 != ""){M=$1;print M"\t"$4}else{print M"\t"$1}}' $Family.tab.sam | sed 's/,/\t/g' | awk '{print $1"\t"$2}' > $Family.list.txt
     samtools faidx $OUT/03-infors-merge/$Family.TE.fa

awk '{for(i=1;i<=1000;i++){ID="id"i;print ID,$1}}' $OUT/03-infors-merge/$Family.TE.fa.fai | \
       awk 'ARGIND==1{A[$1][$2]=1}ARGIND==2{print $2,$1,A[$1][$2]+0}' $Family.list.txt - | \
awk '{
    row=$1;
    col=$2;
    value=$3;
    matrix[row,col]=value;
    rows[row];
    cols[col];
}
END {
    # 打印列索引
    printf("%s ", " ");
    for (col in cols) {
        printf("%s ", col);
    }
    printf("\n");

    # 打印矩阵
    for (row in rows) {
        printf("%s ", row);
        for (col in cols) {
            printf("%s ", matrix[row,col]+0);
        }
        printf("\n");
    }
}' > $Family.matrix.txt

done < $OUT/01-infors.family.txt




