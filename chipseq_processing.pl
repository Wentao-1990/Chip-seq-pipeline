#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

my ($read1,$read2,$sample,$rawdata);
my ($fqsuffix, $genome);
my ($trimfr1,$trimfr2,$trimmr1,$trimmr2,$trimlog);
my ($tophatmapped,$tophatunmapped);
my $tophatout=0;
my $cleardata1=0;
my $middata1=1;
my $middata2=1;
my $cleardata2=0;
my (@read1,@read2);
sub usage{
	die(
		"
	Usage: chipseq_processing.pl -read1 sample1.read1.fq.gz,sample2.read1.fq.gz -read2 sample1.read2.fq.gz,sample2.read2.fq.gz -sample sample1,sample2 -genome indexed genome with bowtie2-build -gtf gtf annotation file -rawdata 1 (raw data)/0 (clear data, default)
	Command: -read1 read1 file path (Must)
		 -read2 read2 file path (Must)
		 -sample bowtie output suffix (Must)
		 -rawdata raw data(1)/ clean data(0)
		 -fqsuffix fq.gz or fastq.gz (Must)
		 -genome indexed genome (Must)
		 \n"
	 )
 }


GetOptions(
	'read1=s'=>\$read1,
	'read2=s'=>\$read2,
	'sample=s'=>\$sample,
	'rawdata=i'=>\$rawdata,
	'fqsuffix=s'=>\$fqsuffix,
	'genome=s'=>\$genome,
);

&usage unless(defined($read1) && defined($sample) && defined($rawdata) && defined($genome) &&defined($fqsuffix));


if($rawdata){
	print "####start trim raw data####\n";
	($trimfr1,$trimfr2,$trimmr1,$trimmr2,$trimlog)=trimrawread($read1,$read2,$sample);
	fastqc($read1,$read2,'fastqc_raw');
	system("nohup multiqc -o fastqc_rawMulti fastqc_raw 2>err.log&");
	while(1){
		if($cleardata1 && $cleardata2 && $middata1==0 && $middata2==0){
			$tophatmapped=alignment($trimfr1,$trimfr2,$sample);
			my $trimfr1_scale=join(',',@$trimfr1);
			my $trimfr2_scale=join(',',@$trimfr2);
			fastqc($trimfr1_scale,$trimfr2_scale,'fastqc_trimmed');
			system("nohup multiqc -o fastqc_trimMulti fastqc_trimmed 2>err.log&");
			last;
		}else{
			$cleardata1=checkfile(@$trimfr1);
			$cleardata2=checkfile(@$trimfr2);
			$middata1=checkfile(@$trimmr1);
			$middata2=checkfile(@$trimmr2);
			#print "$cleardata1\t$cleardata2\t$middata1\t$middata2\n";
			next;
		}
	}
}else{
	@read1=split(/,/,$read1);
	@read2=split(/,/,$read2);
	$tophatmapped=alignment(\@read1,\@read2,$sample);
}





print "####check whether the aglignment is done####\n";

#while(1){
##	if($tophatout){
#		print "Featurecount\n";
#		system("featureCounts -a $gtf  -t 'gene' -F 'GTF' -T 20 -s 0  -B  -C -o allsample_geneid.featureCount @$tophatmapped");
#		last;
#	}else{
#		sleep(100);
#		print "wating for tophat output\n";
#		$tophatout=checkfile(@$tophatunmapped);
#		next;
#	}
#}

sub fastqc{
	my ($r1,$r2,$output)=@_;
	my (@read1,@read2);
	@read1=split(/,/,$r1);
	@read2=split(/,/,$r2);
	system("mkdir $output");
	system("fastqc -o $output -t 2 @read1 @read2");
}
	

sub checkfile{
	my (@filelist)=@_;
	my $tag=0;
	foreach(@filelist){
	        if(-e $_){
			$tag=1;
			next;
		}else{
			$tag=0;
			last;
		}
	}
	return $tag;
}

sub checkfilesize{
	my (@logfile)=@_;
	my $tag=0;

}

sub trimrawread{
	my ($r1,$r2,$sample)=@_;
	my $trim_log;
	my (@trimfr1,@trimfr2,@trimlog);
	my (@trimmr1,@trimmr2);
	my (@read1,@read2,@sample);
	my ($r1_finalname,$r2_finalname);
	my ($r1_midname,$r2_midname);
	my $i;
	@read1=split(/,/,$r1);
	@read2=split(/,/,$r2);
	@sample=split(/,/,$sample);
	if($#read1==$#sample){
		for($i=0;$i<=$#read1;$i++){
			$r1_finalname=basename($read1[$i],".$fqsuffix").'_val_1.fq.gz';
			$r2_finalname=basename($read2[$i],".$fqsuffix").'_val_2.fq.gz';
			$r1_midname=basename($read1[$i],".$fqsuffix").'_trimmed.fq.gz';
			$r2_midname=basename($read2[$i],".$fqsuffix").'_trimmed.fq.gz';
			#print "$r1_finalname\t$r2_finalname\t$r1_midname\t$r2_midname\n";
			$trim_log=basename($read2[$i]).'_trimming_report.txt';
			system "nohup trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 -o  $sample[$i]_trimout -paired $read1[$i] $read2[$i] 2>err.log&\n";
			push @trimfr1,"./$sample[$i]_trimout/$r1_finalname";
			push @trimfr2,"./$sample[$i]_trimout/$r2_finalname";
			push @trimmr1,"./$sample[$i]_trimout/$r1_midname";
			push @trimmr2,"./$sample[$i]_trimout/$r2_midname";
			push @trimlog,"./$sample[$i]_trimout/$trim_log";

		}
	}
	return (\@trimfr1,\@trimfr2,\@trimmr1,\@trimmr2,\@trimlog);
}

sub alignment{
	my ($r1,$r2,$sample)=@_;
	my $i;
	my @mapped;
	my @unmapped;
	my @sample;
	@sample=split(/,/,$sample);
	if($#$r1 == $#$r2){
		#system("mkdir bowtie_alignment");
		for($i=0;$i<=$#$r1;$i++){
			system "mkdir $sample[$i]";
			system "bowtie2 -p 2  -x  $genome -1  $$r1[$i] -2  $$r2[$i] -S $sample[$i]/$sample[$i].bowtiealignment.sam 2>$sample[$i]/$sample[$i].err.log ";
			sleep(10);
			push @mapped,"./bowtie_alignment/$sample[$i]_bowtiealignment.sam";
			#push @unmapped,"./$sample[$i]_tophatout/unmapped.bam";
		}
	}
	return (\@mapped);
}

