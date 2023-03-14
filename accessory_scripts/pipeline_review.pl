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
	"case|c=s"    => \$opts{"case"},
	"pipeline|l=s"=> \$opts{"pipeline"},
	"fpr|f:s"     => \$opts{"fpr"},
	"out|o:s"     => \$opts{"out"},
	"stale|x"     => \$opts{"stale"},
);


validate_options(\%opts);

my $project=$opts{project};
my $case=$opts{case};

my %WFS=(
	WGTS=>{
		bamMergePreprocessing_by_tumor_group=> "alignments_WG.callready",
		bamMergePreprocessing=>                "alignments_WG.callready",
		mutect2_matched_by_tumor_group=>       "calls.mutations",
		variantEffectPredictor_matched_by_tumor_group=>"calls.mutations",
		varscan_by_tumor_group=>"calls.copynumber",
		sequenza_by_tumor_group=>"calls.copynumber",
		delly_matched_by_tumor_group=>"calls.structuralvariants",
		mavis=>"calls.structuralvariants",
		STAR=>"alignments_WT",   ### this will be modifed based on the run to fulldepth vs discrete,
		star_call_ready=>"alignments_WT.callready",
		star_lane_level=>"alignments_WT.lanelevel",
		starfusion=>"calls.fusions",
		rsem=>"calls.expression",
    	###RUO
    	bamMergePreprocessing=>          "alignments_WG.callready",
    	delly_matched=>                  "calls.structuralvariants",
    	mutect2_matched=>                "calls.mutations",
    	sequenza=>                       "calls.copynumber",
    	variantEffectPredictor_matched=> "calls.mutations",
    	varscan=>                        "calls.copynumber",
	},
	WG=>{
		bamMergePreprocessing_by_tumor_group=>"alignments_WG.callready",
		bamMergePreprocessing=>"alignments_WG.callready",
		mutect2_matched_by_tumor_group=>"calls.mutations",
		variantEffectPredictor_matched_by_tumor_group=>"calls.mutations",
		varscan_by_tumor_group=>"calls.copynumber",
		sequenza_by_tumor_group=>"calls.copynumber",
		delly_matched_by_tumor_group=>"calls.structuralvariants",
		mavis=>"calls.structuralvariants",
    	###RUO
    	bamMergePreprocessing=>          "alignments_WG.callready",
    	delly_matched=>                  "calls.structuralvariants",
    	mutect2_matched=>                "calls.mutations",
    	sequenza=>                       "calls.copynumber",
    	variantEffectPredictor_matched=> "calls.mutations",
    	varscan=>                        "calls.copynumber",
	},
	WT=>{
		STAR=>"alignments_WT",   ### this will be modifed based on the run to fulldepth vs discrete,
		star_call_ready=>"alignments_WT.callready",
		starfusion=>"calls.fusions",
		rsem=>"calls.expression",
		mavis=>"calls.structuralvariants",
	}
);





##########################
####. LOAD FPR DATA from this projects
print STDERR "loading data from provenance for project $opts{project}, case $opts{case} for $opts{pipeline} pipeline\n";
my %DATA=loadFPR(%opts);

print STDERR "generating analysis blocks\n";
my %BLOCKS;
%{$BLOCKS{WG}}=pipeline_blocks("WG",%{$DATA{WORKFLOWS}});
%{$BLOCKS{WT}}=pipeline_blocks("WT",%{$DATA{WORKFLOWS}});


print STDERR "generating report\n";
my $report_file=$opts{out} . "/" . $opts{case} . "_" . $opts{pipeline} . "_PIPELINE_REPORT.txt";
(open my $RPT,">",$report_file);
print $RPT "wf_run_id\tn_lanes\tdate\tworkfow\tsamples\n";

my %dot;
for my $PIPE(qw/WG WT/){

	for my $block(keys %{$BLOCKS{$PIPE}}){
		#print "block=$block\n";
		my %BLOCK=%{$BLOCKS{$PIPE}{$block}};
		#print "BLOCK\n", Dumper(%BLOCK) , "\n";

		for my $wfrun(sort{ 
				#$BLOCK{$a} cmp $BLOCK{$b} 
				$DATA{WORKFLOWS}{$a}{date} cmp $DATA{WORKFLOWS}{$b}{date}
			} keys %BLOCK){
			
			my $wf=$DATA{WORKFLOWS}{$wfrun}{wf} || "noworkflow";
			
			my $pad=" " x (50-length($wf));
			my $paddedwf = $wf . $pad;
			
			my @sids=keys %{$DATA{WORKFLOWS}{$wfrun}{sids}};
			my $sids=join(" | ",@sids);
		
			my @limskeys=sort keys %{$DATA{WORKFLOWS}{$wfrun}{limskeys}};
			my $limskeys=join("|",@limskeys);
			my $lanes=scalar @limskeys;
		
			my $date=$DATA{WORKFLOWS}{$wfrun}{date} || "nodate";
			$date=substr($date,0,10);
		
			print $RPT "$wfrun\t$lanes\t$date\t$wf\t$sids\n";

            #### add to the dot file
	    	for my $p_wfrun(keys %{$DATA{WORKFLOWS}{$wfrun}{parents}}){
				## skip if not in the block
				next unless($BLOCK{$p_wfrun});

	    		my $p_wf=$DATA{WORKFLOWS}{$p_wfrun}{wf} || "noworkflow";
				my $pid=$p_wf . "_" . substr($p_wfrun,0,5);
			
				my $cid=$wf . "_" . substr($wfrun,0,5);
				my $pc="$pid -> $cid";
				$dot{$PIPE}{$block}{$pc}++;
	    	}
	    	for my $c_wfrun(keys %{$DATA{WORKFLOWS}{$wfrun}{children}}){
				## skip if not in the block
				next unless($BLOCK{$c_wfrun});

	    		my $c_wf=$DATA{WORKFLOWS}{$c_wfrun}{wf} || "noworkflow";
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
	if(! $opts{case}){
		usage("ERROR : no case provided");
	}
	if(! $opts{fpr}){
		$opts{fpr}="/.mounts/labs/seqprodbio/private/backups/seqware_files_report_latest.tsv.gz";
	}
	if(! $opts{pipeline}){
		usage("ERROR : no pipeline provided");
	}else{
		if($opts{pipeline} ne "WGTS"){
			usage("ERROR: $opts{pipeline} is not a valid pipeline");
		}
	}
	if($opts{out}){
		usage("ERROR: output directory $opts{out} does not exist") if(! -d $opts{out});
	}else{
		$opts{out}=".";
	}
}

sub usage{
	
	print "\npipeline_review.pl : pull analysis blocks from FPR bassed on project/case and predefined pipelines\n";
	print "\nperl pipeline_review.pl [options]\n";
	print "Options are as follows:\n";
	print "\t--project.  The name of the project [required]\n";
	print "\t--case.  The Case (PROJ_nnnn).  expected to be in the File Provenacne Report.  [required]\n";
	print "\t--pipeline. The predefined pipeline for which to generate a report. The current list of pipelines:\n";
	print "\t  WGTS : Whole Genome + Transcriptome pipeline\n";
	print "\t--fpr.  The file provenence report to use, defaults to the system current fpr\n";
	print "\t--out. Output location to place the case report, defaults to the current directory\n";
	print "\t--stale. Allow stale records\n";
	print "\t--help displays this usage message.\n";
	die "\n@_\n\n";
}



### to move to module


sub loadFPR{
	my %opts=@_;
	my $fpr     =$opts{fpr};
	my $project =$opts{project};
	my $case    =$opts{case};


	##################
	#### open and parse the FPR
	########################
	my %WORKFLOWS;my %FILES;  ### these will be returned as sections of a larger hash structure
	(open my $FPR,"gunzip -c $fpr |") || die "unable to open $fpr";
	while(my $rec=<$FPR>){
		my @f=split /\t/,$rec;unshift(@f,"");### add a filed to the begnning to shift indices to 1 based
		my $wf=$f[31];
		#print "$id $wf\n";
	
		#### should this be restricted?? or take all workflow and restrict when they are reported on?
		next unless($WFS{WGTS}{$wf});
		my $platform=$f[23];
		my $platformid=$platform=~/NovaSeq/ ? "NovaSeq" :
			$platform=~/NextSeq/ ? "NextSeq" :
			$platform=~/HiSeq/   ? "HiSeq"   :
			$platform=~/MiSeq/   ? "MiSeq" : "Unknown"; 
				   
		#RESTRICT TO NovaSeq
		next unless($platformid eq "NovaSeq");			   
	
		next unless($f[2] eq $project);
		next unless($f[8] eq $case);
		
		my $wfv=$f[32];
		my $wfrun=$f[37];
		my $inputs=$f[39];     $inputs=~s|vidarr:.*?/file/||g;
		$WORKFLOWS{$wfrun}{wf}=$wf;
		$WORKFLOWS{$wfrun}{wfv}=$wfv;
		$WORKFLOWS{$wfrun}{inputs}=$inputs;
		$WORKFLOWS{$wfrun}{date}=$f[1];
	
	
		my $fid=$f[45];  $fid=~s|vidarr:.*?/file/||g;
		my $file=$f[47];       
		my $md5sum=$f[48];
		$FILES{$fid}={fpath=>$file,md5sum=>$md5sum,wfrun=>$wfrun};
	
		my %info;
		map{
			my($key,$val)=split /=/,$_;
			$info{$key}=$val;
		}split /;/,$f[18];
		my $sid=$case . "_" . $info{geo_tissue_origin} . "_" . $info{geo_tissue_type} . "_" . $info{geo_library_source_template_type};
		if(my $groupid=$info{geo_group_id}){
			$sid.="_${groupid}"
		}
		$WORKFLOWS{$wfrun}{sids}{$sid}={
			tissue_origin=>$info{geo_tissue_origin},
			tissue_type=>$info{geo_tissue_type},
			library_type_=>$info{geo_library_source_template_type},
			group_id=>$info{geo_group_id} || "",
		};
		$WORKFLOWS{$wfrun}{sids}{$sid}{libraries}{$f[14]}++;
		$WORKFLOWS{$wfrun}{sids}{$sid}{limskeys}{$f[57]}++;
		$WORKFLOWS{$wfrun}{limskeys}{$f[57]}++;
	}
	close $FPR;

	#print Dumper(%WORKFLOWS);<STDIN>;
	#print Dumper(%FILES);<STDIN>;

    ########################
	######## identify parent child relationshisp
	###########################
	for my $wfrun(sort keys %WORKFLOWS){
		###########
		####identify workflows with inputs
		#### if there is an input, then use this to define parent child relatioships in with the workflow run information
		#### also use the parent->child info to define an analysis block
		#### all workflow runs that are somehow connected between parents and children should end up in the same block
		if($WORKFLOWS{$wfrun}{inputs}){
			my @inputs=split /,/,$WORKFLOWS{$wfrun}{inputs};  ### adjusted separator from ; to ,
	
			##############
			###### the inputs are file identifiers, which can be used to identify the workflow run that produced the files, ie. the parent workflow
			###########
			for my $fid(@inputs){
				### this is a check for missing parent files
				if(! $FILES{$fid}){
					############# an input exists, but the file was never pulled from file provenance
					$WORKFLOWS{$wfrun}{parents_missing}{$fid}++;
				}else{
					##########. get the parent workflow run
					##########. set as a parent of the current workflow
					##########. set thie curren workflow as a child for the parent workflow
					my $parent_wfrun=$FILES{$fid}{wfrun}; 
					$WORKFLOWS{$wfrun}{parents}{$parent_wfrun}  =$WORKFLOWS{$parent_wfrun}{wf} || "unknown";
	        		$WORKFLOWS{$parent_wfrun}{children}{$wfrun} =$WORKFLOWS{$wfrun}{wf}        || "unknown";
				}
			}
		}
	}
	my %FPR=(WORKFLOWS=>\%WORKFLOWS,FILES=>\%FILES);
	return %FPR;
}

	
sub pipeline_blocks{
	my ($pipeline,%WORKFLOWS)=@_;
	
	#print Dumper(keys %WORKFLOWS);
	
	my %BLOCKS;
	my $blockN=0;
	#######################
	## pass through the set of workflows extracted from FPR
	############	
	for my $A(sort keys %WORKFLOWS){
		###########
		####.use child and parent relationships to define blocks
		#### all workflow runs that are somehow connected between parents and children should end up in the same block
		
		### limit to ONLY WG Workflows
		my $wf=$WORKFLOWS{$A}{wf};
		next unless($WFS{$pipeline}{$wf});

		#print "$wf $wfrun\n";<STDIN>;
		#print Dumper($WORKFLOWS{$A});<STDIN>;
		
		my @p_wfruns=keys %{$WORKFLOWS{$A}{parents}};
		my @c_wfruns=keys %{$WORKFLOWS{$A}{children}};
		my @connected_wfruns=(@p_wfruns,@c_wfruns);
		#print Dumper(@connected_wfruns);<STDIN>;
		#print "parents\n";
		#print Dumper(@p_wfruns);<STDIN>;
		
		
		if(! scalar @connected_wfruns){
			### the workflow had no connections
			### logic : if not in a block already (becasue it was an input to another workflow), then it should be added to a block
			$blockN++;
			$BLOCKS{$A}=$blockN;
		}
		
		
		for my $B(@connected_wfruns){
			#print "A:$A B:$B\n";
			my $wf=$WORKFLOWS{$B}{wf};
			next unless($WFS{$pipeline}{$wf});
			
			
			#####################
			### BLOCK ASSIGNMENT
			### the logic here is expected to produce distinct blocks with all connected workflows
			### the block assigments will change as new workflows are brought in, new blocks are created, or blocks are merged to a new block
			### blocks are incremental values
			### the block structure is WorkflowRun -> BlockNumber

					
			#####
			# logic : if A is in a block, but B not, add B to the A block
			if($BLOCKS{$A} && !$BLOCKS{$B}){
				$BLOCKS{$B}=$BLOCKS{$A};
			}
			# logic : if the B is in a block, but A is not, add A to the B block
			elsif(!$BLOCKS{$A} && $BLOCKS{$B}){
				$BLOCKS{$A}=$BLOCKS{$B};
			}
			# logic : neither A nor B are in a block, create a new block and add both to this
			elsif(!$BLOCKS{$A} && !$BLOCKS{$B}){
				## add both to a new group
				$blockN++;
				$BLOCKS{$A}=$blockN;
				$BLOCKS{$B}=$blockN;
			}
			# logic : both A and B are in a block, but not the same block
			#		  create a new block, and move all from the earlier blocks to the new block
			elsif($BLOCKS{$A} && $BLOCKS{$B} && ($BLOCKS{$A} != $BLOCKS{$B})){
				my $ABlock=$BLOCKS{$A};
				my $BBlock=$BLOCKS{$B};
				$blockN++;  ## the new block
				#print Dumper(%BLOCKS) . "\n";
				map{
					#print "MERGING A:$ABlock B:$BBlock CHECK:$BLOCKS{$_}\n";
					$BLOCKS{$_}=$blockN if( ($BLOCKS{$_}==$ABlock) || ($BLOCKS{$_}==$BBlock));
			  	}keys %BLOCKS;
			}
			# logic : if not the other conditions, they should be in the same group, but lets add a check on this
			else{
				if($BLOCKS{$A} && $BLOCKS{$B} && ($BLOCKS{$A} == $BLOCKS{$B}) ){
					## expect, these are alredy in the same block okay
				}else{
					my $ABlock=$BLOCKS{$A} || "NOBLOCK";
					my $BBlock=$BLOCKS{$B} || "NOBLOCK";
					print STDERR "WARNING : $A and $B BLOCK assignment issue.  ABLOCK=$ABlock, BBLOCK=$BBlock\n";
				}
			}
		}
		
		#print Dumper(%BLOCKS);<STDIN>;
		
		
	}

	#print Dumper(%BLOCKS);<STDIN>;

	################
	### convert the BLOCKS int SETS
	### BLOCKS : workflow run -> block
	### SETS : block -> list of workflow runs
	### get distinct groups and all workflows
	my %SETS;
	map{
		my $set=$BLOCKS{$_};
		$SETS{$set}{$_}=$WORKFLOWS{$_}{wf};
	}keys %BLOCKS;
	
	#print Dumper(%SETS);<STDIN>;
	return %SETS;
}
	


