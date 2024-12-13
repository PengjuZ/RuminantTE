#! /usr/bin/perl
$x=$ARGV[0];
$y=$ARGV[1];
$z=$ARGV[2];
	$LINE=0;
	open (INDEL,"$x");
	while(<INDEL>){
	chomp;
	@A =split (/\s+/ ,$_);
	$LINE+=1;
	}
	$LINE2=0;	
	open (INDEL,"$x");
	while(<INDEL>){
	chomp;	
	@A =split (/\s+/ ,$_);
	$LINE2+=1;
	open (OUT," >> $x.log");
	print OUT "$LINE2/$LINE\n";
	open (OUT," > $x.ID");           		
	print OUT "$A[0]\n";
	open (OUT," > $x.bed");           		
	print OUT "$A[1]	$A[4]	$A[5]\n";	
	system"bedtools getfasta -s -fi $y -bed $x.bed -fo $x.bed.fa";	
	system"seqtk subseq $z $x.ID > $x.ID.fa";
    #system"exonerate -q $x.ID.fa -t $x.bed.fa --model est2genome --querytype dna --targettype dna --showvulgar no --softmaskquery no --softmasktarget yes --showalignment no --showtargetgff yes --showcigar no --score 100 --bestn 1 --percent 60 --refine region";
    system"exonerate -q $x.ID.fa -t $x.bed.fa --model est2genome --querytype dna --targettype dna --showvulgar no --softmaskquery no --softmasktarget yes --showalignment no --showtargetgff yes --showcigar no --score 30 --bestn 1 --percent 50";
	}


