package gsipipeline;
use strict;
use warnings;
use Exporter;


use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION	=	1.00;
@ISA		=	qw/ Exporter /;
@EXPORT		=	qw/ check_file_release get_checks /;




####
#### trying to wrap these up a separate functions calls
my %checks = (
	fastq_1        =>\&fastq_1,
	WG_align_1     =>\&WG_align_1,
	WG_callready_1 =>\&	WG_callready_1,
	WT_align_1     =>\&WT_align_1, 
	WT_callready_1 =>\&WT_callready_1, 
	mutect2_1      =>\&no_release, 
	vep_1          =>\&vep_1, 
	delly_1        =>\&delly_1, 
	mavis_1        =>\&mavis_1,
	varscan_1      =>\&varscan_1,
	sequenza_1     =>\&sequenza_1,
	rsem_1         =>\&rsem_1, 
	starfusion_1   =>\&starfusion_1,
);



sub get_checks { return %checks;  }

sub fastq_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	$status="RELEASE" if($fn=~/\.fastq\.gz$/);
	return ($status,$revised_fn);
}
sub WG_align_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	$status="RELEASE" if($fn=~/\.ba[mi]$/);
	return ($status,$revised_fn);
}
sub WG_callready_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	$status="RELEASE" if($fn=~/\.ba[mi]$/);
	return ($status,$revised_fn);
}
sub WT_align_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	$status="RELEASE" if($fn=~/\.ba[mi]$/);
	return ($status,$revised_fn);
}
sub WT_callready_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	$status="RELEASE" if($fn=~/\.ba[mi]$/);
	return ($status,$revised_fn);
}
	
sub vep_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	$status="RELEASE" if($fn=~/vep/);
	$status="RELEASE" if($fn=~/maf/);
	return ($status,$revised_fn);
}

sub delly_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	### RELEASE the one file
	$status="RELEASE" if($fn=~/delly/ && $fn=~/somatic_filtered/);
	return ($status,$revised_fn);
}

sub mavis_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	### RELEASE ALL
	$status="RELEASE";   	
	return ($status,$revised_fn);
}		

		
sub varscan_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	### release the vcf files
	### release the copy number file
	$status="RELEASE" if($fn=~/vcf$/);
	$status="RELEASE" if($fn=~/copynumber/);
	return ($status,$revised_fn);
}

sub sequenza_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	### release the zipped copy number data from sequnza, with all solutions
	### rename the zipp file to indicate sequenza results
	### modify _somatic -> .somatic
	$status="RELEASE" if($fn=~/results\.zip$/);
	return ($status,$revised_fn);
}
sub rsem_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	### release on the file with results at the end, no modificaitons
	$status="RELEASE" if($fn=~/results$/);
	return ($status,$revised_fn);
}
sub starfusion_1{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	### files need to be renamed, as there is currently NO identifier
	$revised_fn=$sid . "." . $fn;
	### tsv, abridged.tsv, abridged.coding_effect.tsv	
	$revised_fn=$sid . "." . $fn;
	return ($status,$revised_fn);
}		

sub no_release{
	my($fn,$sid)=@_;
	my $status="HOLD";my $revised_fn=$fn;
	### NO DATA RELEASED FROM THIS WORKFLOW
	return ($status,$revised_fn);
}	



sub get_checks { return %checks;  }


############################################
#### this function will take in a registered workflow and a fn and determine if that file shoudl be released
#### there is also the option to generate an alias for the file
########################################




sub check_file_release{
	my($wfset,$fn,$sid)=@_;

	my $status="HOLD";my $revised_fn=$fn;
	### sid is passed for file renaming purposes, if needed
	
	if($wfset eq "fastq_1"){
		###. RELEASE ALL fastq.gz####
		$status="RELEASE" if($fn=~/\.fastq\.gz$/);
	}elsif($wfset eq "WG_align_1" || $wfset eq "WG_callready_1"){
		### RELEASE ONLY bam file and bai index
		$status="RELEASE" if($fn=~/\.ba[mi]$/);

	}elsif($wfset eq "WT_align_1" || $wfset eq "WT_callready_1"){
		### RELEASE ONLY bam file and bai index
		$status="RELEASE" if($fn=~/\.ba[mi]$/);
		
		
	}elsif($wfset eq "mutect2_1"){
		### NO DATA RELEASED FROM THIS WORKFLOW
		
	}elsif($wfset eq "vep_1"){
		### RELEASE vcf file and index with vep in the file name
		### RELEASE zipped maf file
		$status="RELEASE" if($fn=~/vep/);
		$status="RELEASE" if($fn=~/maf/);
		
		### clean up file names
		#($revised_fn=$fn)=~s/filter\.deduped\.realigned\.recalibrated\.mutect2/mutect2/;
		
	}elsif($wfset eq "delly_1"){
		### RELEASE the one file
		$status="RELEASE" if($fn=~/delly/ && $fn=~/somatic_filtered/);
		
		### get rid of the first _somatic which comes from the olive
		#($revised_fn=$fn)=~s/_+somatic//;
		## there are cases where files are named with WG_.somatic.filtered.delly.merged.vcf, correcting to get rid of that. underscore
		#$revised_fn=~s/WG_\./WG\./;
		#$revised_fn=~s/somatic_filtered\.delly\.merged\.vcf/delly\.somatic\.vcf/;
	
	
	}elsif($wfset eq "mavis_1"){
		### RELEASE ALL
		$status="RELEASE";   
		
	}elsif($wfset eq "varscan_1"){
		### release the vcf files
		### release the copy number file
		$status="RELEASE" if($fn=~/vcf$/);
		$status="RELEASE" if($fn=~/copynumber/);
		
	}elsif($wfset eq "sequenza_1"){
		### release the zipped copy number data from sequnza, with all solutions
		### rename the zipp file to indicate sequenza results
		### modify _somatic -> .somatic
		$status="RELEASE" if($fn=~/results\.zip$/);
	
	    #print "filename = $fn\n";
		
		#$revised_fn=$fn;
		#$revised_fn=~s/_somatic/\.varscan\.somatic/;
		#$revised_fn=~s/_results/.sequenza\.results/;
		
		#($revised_fn=$fn)=~s/_somatic/\.somatic/;
		#($revised_fn=$fn)=~s/_results/.sequenza_results/;
		
		#print "revised   = $revised_fn\n";
		
	}elsif($wfset eq "rsem_1"){
		### release on the file with results at the end, no modificaitons
		$status="RELEASE" if($fn=~/results$/);
		#$revised_fn=$fn;
		#$revised_fn=~s/genes/RSEM\.genes/;
		#$revised_fn=~s/isoforms/RSEM\.isoforms/;
		
	}elsif($wfset eq "starfusion_1"){
		### release only the tsv files.  All of them?
		$status="RELEASE" if($fn=~/tsv$/);
		### files need to be renamed, as there is currently NO identifier
		$revised_fn=$sid . "." . $fn;
		### tsv, abridged.tsv, abridged.coding_effect.tsv
		
	}else{
		print STDERR "$wfset NOT REGISTERED\n";
	}
	
	return ($status,$revised_fn);
}
