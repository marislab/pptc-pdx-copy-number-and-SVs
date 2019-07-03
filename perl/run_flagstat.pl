use Parallel::ForkManager;
use Data::Dumper;
use File::Slurp;

my $file_location = '/mnt/isilon/maris_lab/target_nbl_ngs/rathik/ATRX_exon_coverage/flagstat_output.txt';
open my $file, '>>', $file_location or die $!;

my @samples = read_file '/mnt/isilon/maris_lab/target_nbl_ngs/rathik/ATRX_exon_coverage/samples.txt';
my $bamdir = '/mnt/isilon/maris_lab/target_nbl_ngs/PPTC-PDX-genomics/bcm-processed/wes-bam';
my $outdir = '/mnt/isilon/maris_lab/target_nbl_ngs/rathik/ATRX_exon_coverage/flagstat_results';

# start time
my $start_run = time();

# run parallel jobs
my $pm=new Parallel::ForkManager(10);

foreach (@samples)
{	
	$pm->start and next;

		chomp($_);
		$sample = $_;
		$sname = $sample;
		$sname =~ s/-human.bam|.hybrid_hg19.realigned.recal.bam/_coverage.txt/;

		print $sample,"\n";
		print $sname,"\n";
				
		# my $bedcov = "samtools flagstat $bamdir/$sample > $outdir/$sname";
		# print $bedcov,"\n";
		# system($bedcov);

		$fname = $sample;
		$fname =~ s/-human.bam|.hybrid_hg19.realigned.recal.bam//;
		my $flagstat = `samtools flagstat $bamdir/$sample | head -n1 | sed 's/ .*//g'`;
		# print $flagstat,"\n";
		# $var = system($flagstat);
		$res = "$fname\t$flagstat";
		print $file "$res\n";

	$pm->finish;
}

$pm->wait_all_children;

# end time
my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";
