package gsipipeline;
use strict;
use warnings;
use Exporter;
use JSON::PP;
use Data::Dumper;


use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION	=	1.00;
@ISA		=	qw/ Exporter /;
@EXPORT		=	qw/ check_file_release get_checks loadFPR pipeline_blocks/;




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




############################################
### these functions are for pipeline review
############################

sub loadFPR{
	my %opts=@_;
	my $fpr     =$opts{fpr};
	my $project =$opts{project};
	#my $case    =$opts{case};
	my %caselist = map{$_=>1} @{$opts{caselist}};  ### plase list of cases into a hash, for easy checking
	#print Dumper(%caselist);exit;
	my $pipeline = $opts{pipeline};
	my $PIPELINES= $opts{pipelines};

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
		next unless($$PIPELINES{$pipeline}{Workflows}{$wf});
		my $platform=$f[23];
		my $platformid=$platform=~/NovaSeq/ ? "NovaSeq" :
			$platform=~/NextSeq/ ? "NextSeq" :
			$platform=~/HiSeq/   ? "HiSeq"   :
			$platform=~/MiSeq/   ? "MiSeq" : "Unknown"; 
				   
		#RESTRICT TO NovaSeq
		next unless($platformid eq "NovaSeq");			   
	
		next unless($f[2] eq $project);
		my $case=$f[8];
		next unless($caselist{$case});
		
		my $wfv=$f[32];
		my $wfrun=$f[37];
		my $inputs=$f[39];     $inputs=~s|vidarr:.*?/file/||g;
		$WORKFLOWS{$case}{$wfrun}{wf}=$wf;
		$WORKFLOWS{$case}{$wfrun}{wfv}=$wfv;
		$WORKFLOWS{$case}{$wfrun}{inputs}=$inputs;
		$WORKFLOWS{$case}{$wfrun}{date}=$f[1];
	
	
		my $fid=$f[45];  $fid=~s|vidarr:.*?/file/||g;
		my $file=$f[47];       
		my $md5sum=$f[48];
		$FILES{$case}{$fid}={fpath=>$file,md5sum=>$md5sum,wfrun=>$wfrun};
	
		my %info;
		map{
			my($key,$val)=split /=/,$_;
			$info{$key}=$val;
		}split /;/,$f[18];
		my $sid=$case . "_" . $info{geo_tissue_origin} . "_" . $info{geo_tissue_type} . "_" . $info{geo_library_source_template_type};
		if(my $groupid=$info{geo_group_id}){
			$sid.="_${groupid}"
		}
		$WORKFLOWS{$case}{$wfrun}{sids}{$sid}={
			tissue_origin=>$info{geo_tissue_origin},
			tissue_type=>$info{geo_tissue_type},
			library_type_=>$info{geo_library_source_template_type},
			group_id=>$info{geo_group_id} || "",
		};
		$WORKFLOWS{$case}{$wfrun}{sids}{$sid}{libraries}{$f[14]}++;
		$WORKFLOWS{$case}{$wfrun}{sids}{$sid}{limskeys}{$f[57]}++;
		$WORKFLOWS{$case}{$wfrun}{limskeys}{$f[57]}++;
	}
	close $FPR;

	#print Dumper(%WORKFLOWS);<STDIN>;
	#print Dumper(%FILES);<STDIN>;

    ########################
	######## identify parent child relationshisp
	###########################
	for my $case(sort keys %WORKFLOWS){
		for my $wfrun(sort keys %{$WORKFLOWS{$case}}){
			if($wfrun eq "14522cc0c2f5fa873af81cbe2dc8afea1c1ba6f7dc2c1955b73abc46dd9fabf2"){
			}
			###########
			####identify workflows with inputs
			#### if there is an input, then use this to define parent child relatioships in with the workflow run information
			if($WORKFLOWS{$case}{$wfrun}{inputs}){
				my @inputs=split /,/,$WORKFLOWS{$case}{$wfrun}{inputs};  ### adjusted separator from ; to ,
				##############
				###### the inputs are file identifiers, which can be used to identify the workflow run that produced the files, ie. the parent workflow
				###########
				for my $fid(@inputs){
					### this is a check for missing parent files
					if(! $FILES{$case}{$fid}){
						############# an input exists, but the file was never pulled from file provenance
						$WORKFLOWS{$case}{$wfrun}{parents_missing}{$fid}++;
					}else{
						##########. get the parent workflow run
						##########. set as a parent of the current workflow
						##########. set thie curren workflow as a child for the parent workflow
						my $parent_wfrun=$FILES{$case}{$fid}{wfrun}; 
						$WORKFLOWS{$case}{$wfrun}{parents}{$parent_wfrun}  =$WORKFLOWS{$case}{$parent_wfrun}{wf} || "unknown";
	        			$WORKFLOWS{$case}{$parent_wfrun}{children}{$wfrun} =$WORKFLOWS{$case}{$wfrun}{wf}        || "unknown";
					}
				}
			}
		}
	}
	my %FPR=(WORKFLOWS=>\%WORKFLOWS,FILES=>\%FILES);
	return %FPR;
}

	
sub pipeline_blocks{
	my ($pipeline,$WORKFLOWS,$PIPELINES)=@_;
	my %WFBLOCKS;
	my $blockN=0;
	#######################
	## pass through the set of workflows extracted from FPR
	############	
	for my $A(sort keys %$WORKFLOWS){
		###########
		####.use child and parent relationships to define blocks
		#### all workflow runs that are somehow connected between parents and children should end up in the same block
		
		### limit to ONLY WG Workflows
		my $wf=$$WORKFLOWS{$A}{wf};
		next unless($$PIPELINES{$pipeline}{Workflows}{$wf});

		#print "$wf $wfrun\n";<STDIN>;
		
		
		my @p_wfruns=keys %{$$WORKFLOWS{$A}{parents}};
		my @c_wfruns=keys %{$$WORKFLOWS{$A}{children}};
		my @connected_wfruns=(@p_wfruns,@c_wfruns);
		#print Dumper(@connected_wfruns);<STDIN>;
		#print "parents\n";
		#print Dumper(@p_wfruns);<STDIN>;
		
		
		if(! scalar @connected_wfruns){
			### the workflow had no connections
			### logic : if not in a block already (becasue it was an input to another workflow), then it should be added to a block
			print STDERR "registering unconnected run\n";
			print "run=$A\n" . Dumper($$WORKFLOWS{$A});<STDIN>;
			$blockN++;
			$WFBLOCKS{$A}=$blockN;
		}
		
		
		for my $B(@connected_wfruns){
			#print "A:$A B:$B\n";
			my $wf=$$WORKFLOWS{$B}{wf};
			next unless($$PIPELINES{$pipeline}{Workflows}{$wf});
			
			
			#####################
			### BLOCK ASSIGNMENT
			### the logic here is expected to produce distinct blocks with all connected workflows
			### the block assigments will change as new workflows are brought in, new blocks are created, or blocks are merged to a new block
			### blocks are incremental values
			### the block structure is WorkflowRun -> BlockNumber

					
			#####
			# logic : if A is in a block, but B not, add B to the A block
			if($WFBLOCKS{$A} && !$WFBLOCKS{$B}){
				$WFBLOCKS{$B}=$WFBLOCKS{$A};
			}
			# logic : if the B is in a block, but A is not, add A to the B block
			elsif(!$WFBLOCKS{$A} && $WFBLOCKS{$B}){
				$WFBLOCKS{$A}=$WFBLOCKS{$B};
			}
			# logic : neither A nor B are in a block, create a new block and add both to this
			elsif(!$WFBLOCKS{$A} && !$WFBLOCKS{$B}){
				## add both to a new group
				$blockN++;
				$WFBLOCKS{$A}=$blockN;
				$WFBLOCKS{$B}=$blockN;
			}
			# logic : both A and B are in a block, but not the same block
			#		  create a new block, and move all from the earlier blocks to the new block
			elsif($WFBLOCKS{$A} && $WFBLOCKS{$B} && ($WFBLOCKS{$A} != $WFBLOCKS{$B})){
				my $ABlock=$WFBLOCKS{$A};
				my $BBlock=$WFBLOCKS{$B};
				$blockN++;  ## the new block
				#print Dumper(%WFBLOCKS) . "\n";
				map{
					#print "MERGING A:$ABlock B:$BBlock CHECK:$WFBLOCKS{$_}\n";
					$WFBLOCKS{$_}=$blockN if( ($WFBLOCKS{$_}==$ABlock) || ($WFBLOCKS{$_}==$BBlock));
			  	}keys %WFBLOCKS;
			}
			# logic : if not the other conditions, they should be in the same group, but lets add a check on this
			else{
				if($WFBLOCKS{$A} && $WFBLOCKS{$B} && ($WFBLOCKS{$A} == $WFBLOCKS{$B}) ){
					## expect, these are alredy in the same block okay
				}else{
					my $ABlock=$WFBLOCKS{$A} || "NOBLOCK";
					my $BBlock=$WFBLOCKS{$B} || "NOBLOCK";
					print STDERR "WARNING : $A and $B BLOCK assignment issue.  ABLOCK=$ABlock, BBLOCK=$BBlock\n";
				}
			}
		}
		
		#print Dumper(%WFBLOCKS);<STDIN>;
		
		
	}

	#print Dumper(%WFBLOCKS);<STDIN>;
    #print Dumper(values %WFBLOCKS);<STDIN>;
	## rename the blocks
	my @blocks = do { my %seen; grep { !$seen{$_}++ } values %WFBLOCKS };
	my $n=0;
	my %newblocks;
	map{
		$n++;
		$newblocks{$_}="${pipeline}.$n";
	}sort @blocks;
	#print Dumper(%newblocks);<STDIN>;
	################
	### convert the WFBLOCKS int BLOCKS
	### WFBLOCKS : workflow run -> block
	### BLOCKS : block -> list of workflow runs
	### get distinct groups and all workflows
	my %BLOCKS;
	map{
		my $block=$WFBLOCKS{$_};
		my $newblock=$newblocks{$block};
		$BLOCKS{$newblock}{$_}=$$WORKFLOWS{$_}{wf};
	}sort keys %WFBLOCKS;
	
	#print Dumper(%BLOCKS);<STDIN>;
	return %BLOCKS;
}
