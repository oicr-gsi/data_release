#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Getopt::Long;


### pull libraries from the same directory at this file
use lib dirname (__FILE__);
use gsipipeline;


my %opts = ();
GetOptions(
	"project|p=s" => \$opts{"project"},
	"fpr|f:s"     => \$opts{"fpr"},
	"out|o:s"     => \$opts{"out"},
	"runfile|r:s" => \$opts{"runfile"},
);
validate_options(\%opts);

my %checks=get_checks();




#################
### workflows+ versions are registered here and assigned to a more general name
### multiple workflows may map to the same general name, if so they should all have the same files, at least those being released
### these key words will be used in the check file release section
### will move to getting this from an external file where new workflows can be registered
################
my %registered_workflow_sets=(
	delly_1           =>[qw/delly_matched_by_tumor_group_2.3.1/],
	fastq_1           =>[qw/bcl2fastq_3.1.3/],
	mavis_1           =>[qw/mavis_3.0.3/],
	mutect2_1         =>[qw/mutect2_matched_by_tumor_group_1.0.4/],
	rsem_1            =>[qw/rsem_1.0.0/],
	sequenza_1        =>[qw/sequenza_by_tumor_group_1.0.0/],
	starfusion_1      =>[qw/starfusion_2.0.1/],
	WG_align_1        =>[qw/bwaMem_1.0.0/],
	WG_callready_1    =>[qw/bamMergePreprocessing_by_tumor_group_2.0.3/],
	WT_align_1        =>[qw/star_lane_level_2.0.4/],
	WT_callready_1    =>[qw/star_call_ready_2.0.4/],
	vep_1             =>[qw/variantEffectPredictor_matched_by_tumor_group_2.1.7 variantEffectPredictor_matched_by_tumor_group_2.1.8/],
	varscan_1         =>[qw/varscan_by_tumor_group_2.2.4/],					
);
my %registered_workflows;
for my $set(keys %registered_workflow_sets){
	for my $wf(@{$registered_workflow_sets{$set}}){
		$registered_workflows{$wf}=$set;
	}
}
my %pipeline_sets=(

	WGTS=>qw/fastq_1 WG_align_1 mutect2_1 varscan_1 delly_1 vep_1 sequenza_1 mavis_1 rsem_1 starfusion_1/,   ## WGTS release excludes lane level alignments
);


#print Dumper(%registered_workflows);
#exit;







#print Dumper(%wfruns);	





#######################
####. MAIN
#######################
my %wfruns;
print STDERR "reading workflow runs from $opts{runfile}\n";
(open my $RUNS,"<","$opts{runfile}") || die "could notopen $opts{runfile}";
while(<$RUNS>){
	chomp;
	$wfruns{$_}{listcount}++;
}


##########################################
###### COLLECTING FILES FROM FPR
##########################################
my %dataset;
my $wfruncount=scalar keys %wfruns;
print STDERR "retrieving data for project $opts{project} from $wfruncount workflow runs\n";
print STDERR "parsing fpr $opts{fpr}\n";
(open my $FPR,"gunzip -c $opts{fpr} |") || die "unable to open fpr at $opts{fpr}";
while(my $rec=<$FPR>){
	#print "$rec";<STDIN>;
	my @f=split /\t/,$rec;unshift(@f,"");
	my $project=$f[2];
	my $case=$f[8];
	#next unless( $cases{$case} && $cases{$case}{requested});
	
	my $wf=$f[31];
	my $wfv=$f[32];
	my $wfrun=$f[37];
	next unless($wfruns{$wfrun});
	$wfruns{$wfrun}{wf}=$wf;
	$wfruns{$wfrun}{wfv}=$wf;
	my $workflow="${wf}_${wfv}";


	
	my $fid=$f[45];
	my $fpath=$f[47];
	my $fn=basename($fpath);
	my $md5sum=$f[48];
	my $limskey=$f[57];
	
	
	my %info;
	map{
		my($key,$val)=split /=/,$_;
		$info{$key}=$val;
	}split /;/,$f[13];
	
	my $sid=$case . "_" . $info{geo_tissue_origin} . "_" . $info{geo_tissue_type} . "_" . $info{geo_library_source_template_type};
	$sid.="_$info{geo_group_id}" if($info{geo_group_id});
	print "$sid\n";
	
	
	
	

	### check to see if this file is released, based on wf + wf ver

	my $regwf=$registered_workflows{$workflow} || 0;
	if($regwf){
		my($status,$release_fn)=$checks{$regwf}->($fn,$sid);
		if($status eq "RELEASE"){
			$dataset{RELEASED}{$fpath}{fn}=$release_fn;
			$dataset{RELEASED}{$fpath}{rpt}="$project\t$case\t$regwf\t$workflow\t$wfrun\t$fn\t$fid\t$release_fn\t$md5sum\t$fpath";
		}else{
			$dataset{UNRELEASED}{$fpath}{rpt}="$project\t$case\t$regwf\t$workflow\t$wfrun\t$fn\t$fid\t$release_fn\t$md5sum\t$fpath";
		}
		
	}else{
		print STDERR "workflow $workflow is not registered, workflow run = $wfrun\n";
	}
	
	#$wfruns{$wfrun}{$fswid}{path}=$f[47];
	#$wfruns{$wfrun}{$fswid}{fn}=basename($f[47]);
	#$wfruns{$wfrun}{$fswid}{md5sum}=$f[48];
}


#### the process will create a project directory in the output space
### if this exists already, it should be removed
my $project_dir="$opts{out}/$opts{project}";
mkdir($project_dir) unless(-d $project_dir);
print STDERR "opening report file for Released files\n";
my $released_file=$opts{out} . "/" . $opts{project} . "_release_report.txt";
(open my $RPT,">",$released_file) || die "unable to open $released_file";
print $RPT "Project\tCase\tpipestep\tworkflow\twfrun\tfn\tfid\trelease_fn\tmd5sum\tpath\n";

for my $fpath(sort keys %{$dataset{RELEASED}}){
	my $fn=$dataset{RELEASED}{$fpath}{fn};
	my $line=$dataset{RELEASED}{$fpath}{rpt};
	
	
	my $ln=$opts{out} . "/" . $fn;
	if(-e $ln){
		print STDERR "ERROR file already linked out $ln";<STDIN>;
	}else{
		symlink($fpath,$ln) || die "coulnd not link file to $ln from $fpath";
	}
	
	print $RPT "$line\n";
}
close $RPT;

print STDERR "opening report file for Unreleased files\n";
my $unreleased_file=$opts{out} . "/" . $opts{project} . "_UNRELEASED_report.txt";
(open my $UNRELEASED,">",$unreleased_file) || die "unable to open $unreleased_file";
print $UNRELEASED "Project\tCase\tpipestep\tworkflow\twfrun\tfn\tfid\trelease_fn\tmd5sum\tpath\n";
for my $fpath(sort keys %{$dataset{UNRELEASED}}){
	my $line=$dataset{UNRELEASED}{$fpath}{rpt};
	print $UNRELEASED "$line\n";
}
close $UNRELEASED;








sub validate_options{
	my ($opts)=@_;
	usage("Help requested.") if($opts{help});
	if(! $opts{project}){
		usage("ERROR : no project provided");
	}
	if(! $opts{cases}){
		$opts{cases}="ALL";
	}else{
		### check for a comma separated list of cases
		if($opts{cases}=~/\s/){
			usage("ERROR : case list should be a comma separated with no whitespace")
		}
	}	
	if(! $opts{runfile}){
		usage("ERROR : no runfile provided");
	}
	if(! $opts{fpr}){
		$opts{fpr}="/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz";
	}
	if($opts{out}){
		usage("ERROR: output directory $opts{out} does not exist") if(! -d $opts{out});
		
	}else{
		$opts{out}=".";
	}
}

sub usage{
	
	print "\npipeline_release.pl : extracts files and locations from file provenance for a set of workflow runs based the pipeline output\n";
	print "\npipeline_release.pl [options]\n";
	print "Options are as follows:\n";
	print "\t--project.  The name of the project [required]\n";
	#print "\t--cases.  A comma separated list of cases, with no whitespace.  [optional]\n";
	print "\t--runfile.  List of workflow runs to link out.  First column must be a valid workflow run accession, expected to be in the provided/default FPR\n";
	print "\t--fpr.  The file provenence report to use, defaults to the system current fpr\n";
	print "\t--out. Output directory. where to put the links.  This will be organized by Project/Case.\n";
	print "\t--help displays this usage message.\n";
	die "\n@_\n\n";
}


