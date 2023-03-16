#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use JSON::PP;


### pull libraries from the same directory at this file
use lib dirname (__FILE__);
use gsipipeline;



my %opts = ();
GetOptions(
	"project|p=s"  => \$opts{"project"},
	"cases|c=s"    => \$opts{"cases"},
	"pipeline|l=s" => \$opts{"pipeline"},
	"fpr|f:s"      => \$opts{"fpr"},
	"outdir|d:s"   => \$opts{"outdir"},
	"outfile|o:s"  => \$opts{"outfile"},
	"stale|x"      => \$opts{"stale"},
	"graph|g=s"    => \$opts{"graph"},
	"config|f=s"   => \$opts{"config"},
	"help|h"       => \$opts{"help"},
);
validate_options(\%opts);

## load the pipeline configuration information
### this gives descriptions of the various pipelines + workflow names that are registered with each
print STDERR "loading the pipeline configuration file\n";
(open my $JSON,"<",$opts{config}) || die "could not load the pipeline configuration file $opts{config}";
#my $json_text=join("",<$JSON>);
my $json_text=`cat $opts{config}`;
close $JSON;


my $PIPELINES=decode_json($json_text);
### add this to opts for each passing around
$opts{pipelines}=$PIPELINES;




my $project=$opts{project};
my $case=$opts{case};  ### will convert to cases, a comma separated list, so all can be collected in a single report
my $pipeline=$opts{pipeline};
### check that the pipeline is supports
if(! $$PIPELINES{$pipeline}){
	my $supported=join("|",keys %$PIPELINES);
	usage("Requested pipeline $pipeline is not supported. Supported pipelines are $supported")
}




##########################
####. LOAD FPR DATA from this projects
print STDERR "loading data from provenance for project $opts{project}, cases $opts{cases} for $opts{pipeline} pipeline\n";
my %DATA=loadFPR(%opts);
#print Dumper(%DATA);exit;



for my $case(@{$opts{caselist}}){
	print STDERR "processing case $case\n";
	
	#### check if there is data for this case
	if(! $DATA{WORKFLOWS}{$case}){
		print STDERR "WARNING : no workflows workflows for this pipeline were found for case $case\n";
		next;
	}
    my %WORKFLOWS=%{$DATA{WORKFLOWS}{$case}};
	#print STDERR Dumper(keys %WORKFLOWS);
	
	print STDERR "generating analysis blocks\n";
	my %BLOCKS;
	#print "BLOCKS=". Dumper(%BLOCKS);<STDIN>;

	if($opts{pipeline} eq "WGTS"){
		### to do, merge to a single set of blocks
		%{$BLOCKS{WG}}=pipeline_blocks("WG",\%WORKFLOWS,$PIPELINES);
		%{$BLOCKS{WT}}=pipeline_blocks("WT",\%WORKFLOWS,$PIPELINES);
	}

	if($opts{pipeline}eq "WG"){
		%{$BLOCKS{WG}}=pipeline_blocks("WG",\%WORKFLOWS,$PIPELINES);
	}

	if($opts{pipeline}eq "WT"){
		%{$BLOCKS{WT}}=pipeline_blocks("WT",\%WORKFLOWS,$PIPELINES);
	}

	print STDERR "generating report\n";
	my $report_file=$opts{outdir} . "/" . $case . "_" . $opts{pipeline} . "_PIPELINE_REPORT.txt";
	(open my $RPT,">",$report_file);
	print $RPT "WORKFLOW RUN REPORT for Project:$project\n";
	print $RPT "Pipeline:$pipeline , $$PIPELINES{$pipeline}{Title}\n$$PIPELINES{$pipeline}{Description}\n\n";

	print $RPT "Case=$case\n";

	my %dot;
	for my $PIPE(qw/WG WT/){
		my $blockcount=scalar keys %{$BLOCKS{$PIPE}};
		print STDERR "$blockcount BLOCKS for $PIPE pipeline\n";
		
		for my $block(sort keys %{$BLOCKS{$PIPE}}){
			#print "block=$block\n";
			my %BLOCK=%{$BLOCKS{$PIPE}{$block}};
			#print "BLOCK\n", Dumper(%BLOCK) , "\n";
			print $RPT "#BLOCK $block\n";
			print $RPT "wf_run_id\tn_lanes\tdate\tworkflow\tsamples\tlimskeys\n";
			for my $wfrun(sort{ 
					#$BLOCK{$a} cmp $BLOCK{$b} 
					$WORKFLOWS{$a}{date} cmp $WORKFLOWS{$b}{date}
				} keys %BLOCK){
				#print Dumper($DATA{WORKFLOWS}{$wfrun});<STDIN>;
				my $wf=$WORKFLOWS{$wfrun}{wf} || "noworkflow";
			
				my $pad=" " x (50-length($wf));
				my $paddedwf = $wf . $pad;
			
				my @sids=sort keys %{$WORKFLOWS{$wfrun}{sids}};
				my $sids=join(" | ",@sids);
		
				my @limskeys=sort keys %{$WORKFLOWS{$wfrun}{limskeys}};
				my $limskeys=join("|",@limskeys);
				my $lanes=scalar @limskeys;
		
				my $date=$WORKFLOWS{$wfrun}{date} || "nodate";
				$date=substr($date,0,10);
		
				print $RPT "$wfrun\t$lanes\t$date\t$paddedwf\t$sids\t$limskeys\n";

            	#### add to the dot file
	    		for my $p_wfrun(keys %{$WORKFLOWS{$wfrun}{parents}}){
					## skip if not in the block
					next unless($BLOCK{$p_wfrun});

	    			my $p_wf=$WORKFLOWS{$p_wfrun}{wf} || "noworkflow";
					my $pid=$p_wf . "_" . substr($p_wfrun,0,5);
			
					my $cid=$wf . "_" . substr($wfrun,0,5);
					my $pc="$pid -> $cid";
					$dot{$PIPE}{$block}{$pc}++;
	    		}
	    		for my $c_wfrun(keys %{$WORKFLOWS{$wfrun}{children}}){
					## skip if not in the block
					next unless($BLOCK{$c_wfrun});

	    			my $c_wf=$WORKFLOWS{$c_wfrun}{wf} || "noworkflow";
					my $cid=$c_wf . "_" . substr($c_wfrun,0,5);
			
					my $pid=$wf . "_" . substr($wfrun,0,5);
					my $pc="$pid -> $cid";
					$dot{$PIPE}{$block}{$pc}++;
	    		}
			}
			print $RPT "\n";
			#print "dot\n" . Dumper($dot{$block});<STDIN>;
		}
	}	
	close $RPT;
}




#print Dumper(%dot);
#for my $block(sort keys %dot){
#	(open my $DOT,">","graph.$block.dot") || die "could not open dot file";
#	print $DOT "digraph {\n";
#	for my $key(sort keys %{$dot{$block}}){
#		print $DOT " $key;\n"
#	}
#	print $DOT "}";
#	#close $DOTT:
#}















sub validate_options{
	my ($opts)=@_;
	usage("Help requested.") if($opts{help});
	if(! $opts{project}){
		usage("ERROR : no project provided");
	}
	if(! $opts{cases}){
		usage("ERROR : no cases provided");
	}else{
		### cases is expected to be comma separated, split this by comma to generate an array of cases
		@{$opts{caselist}}=split /,/,$opts{cases};
	}
	if(! $opts{fpr}){
		$opts{fpr}="/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz";
	}
	if(! $opts{pipeline}){
		usage("ERROR : no pipeline provided");
	}
	if($opts{outdir}){
		usage("ERROR: output directory $opts{outdir} does not exist") if(! -d $opts{outdir});
	}else{
		$opts{outdir}=".";
	}
	if(!$opts{config}){
		my $dir=dirname (__FILE__);
		$opts{config}="$dir/gsiPipelines.json";
	}elsif( ! -e $opts{config} ){
		usage("ERROR: pipeline configuration file $opts{config} not found");
	}
}

sub usage{
	
	print "\npipeline_review.pl : pull analysis blocks from FPR bassed on project/case and predefined pipelines\n";
	print "\nperl pipeline_review.pl [options]\n";
	print "Options are as follows:\n";
	print "\t--project.  The name of the project [required]\n";
	print "\t--cases. A list of cases (PROJ_nnnn), expected to be in the File Provenence Report. This should be comma-seprated with no spaces  [required]\n";
	print "\t--pipeline. The predefined pipeline for which to generate a report. The current list of pipelines:\n";
	print "\t--config. The configuration file (json) that defines the pipelines. Defaults to gsiPipelines.json bundled with this script\n";
	print "\t--fpr.  The file provenence report to use, defaults to the system current fpr\n";
	print "\t--outdir. Output location to place the case report, defaults to the current directory\n";
	print "\t--stale. Allow stale records\n";
	print "\t--help displays this usage message.\n";
	die "\n@_\n\n";
}







	


