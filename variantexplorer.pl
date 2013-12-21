#!/usr/local/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Getopt::Std;
use vars qw/$opt_t/;

getopts('t');
# infile is the first command line argument and outfile the second
my $infile  = $ARGV[0];
my $outfile = $ARGV[1];

# printing the command line arguments
print "Input  file : ", $infile, "\n";
print "Output file : ", $outfile, "\n";
if ($opt_t) { print "Recognized option -t \n" };

# setting up de registry to connect to Ensembl
my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
    );

# adaptor for Core 
my $sa = $reg->get_adaptor("human", "core", "slice");
# adaptor for Variation
my $vfa = $reg->get_adaptor("human", "variation", "variationfeature");
# adaptor for Regulation
my $regfeat_adaptor = $reg->get_adaptor('Human', 'funcgen', 'regulatoryfeature');

## for polyphen and sift and consequences
#my $transcript_adaptor = $reg->get_adaptor('homo_sapiens', 'core', 'transcript'); #get the adaptor to get the Transcript from the database
#my $trv_adaptor = $reg->get_adaptor('homo_sapiens', 'variation', 'transcriptvariation'); #get the adaptor to get TranscriptVariation objects


my @populations;

# read config.txt
open(CONFIGFILE,"config.txt") || die "Could not open configuration file config.txt";
while(<CONFIGFILE>){ # loop through config file line by line 
    chomp;  # remove newline character from end of line
    my $stripped = $_; 
    $stripped =~ s/\s+//g; # remove spaces from line
    if (length($stripped)>1 && ( substr($stripped,0,1) ne "#")) { push(@populations, $stripped); } # add line to populations array
}
close(CONFIGFILE);

foreach my $pop (@populations) {
    print "Using population: ",$pop, "\n"; # print all populations used
}

open(FILE,$infile) || die "Could not open the file ${infile} \n";
my $debug = 0;
my @val;
my @all;
my $i;

open OUTFILE, ">", $outfile or die $!;  ## open outfile in append mode

while(<FILE>){ # loop through each line in the input VCF file
	chomp; # remove newline from the end of the line

	# this part processes the comment lines
	# first it checks whether the line starts with the # symbol
	# in case option -t is given and the line starts with #CHROM, 
	# the line is split based on tabs, and results put into the array @val
	# for each of the data items that we will add to the file
	# a header is added to the $val[7] so that they appear directly after the INFO field
	# all headers are separated by a tab
	# in case option -t is not given, the comment line is simply copied to the output
	if ( /^#/ ) {
		if ($opt_t && /^#CHROM/ ) {
			@val = split(/\t/,$_);
			$val[7] .= "\tDNASE1\tHISTONE\tPOLYMERASE\tTFBS\tDNASE1_CELLTYPES\tHISTONE_CELLTYPES\tPOLYMERASE_CELLTYPES\tTFBS_CELLTYPES\tGENE_ID\tGENE_NAME\tBIOTYPE";
			foreach my $pop (@populations) {  ## add all populations to column 7
				$val[7] = $val[7]."\t${pop}_REF\t${pop}_REFFREQ\t${pop}_ALT\t${pop}_ALTFREQ";
			}
			print OUTFILE join("\t",@val);
			print OUTFILE "\n";
		} else {
			print OUTFILE $_,"\n";
		}
		next; ## go to next while iteration
	}

    # here non-comment lines are processed
	# line split on tabs and stored into the array @val
	# chromosome and position are the 0th and 1st value
	
	@val = split(/\t/,$_);
	my $chromosome=$val[0];
	my $position=$val[1];
	
	# get the slice for this position
	my $slice = $sa->fetch_by_region('chromosome', $chromosome, $position, $position);
    
    ## retrieve gene ID, gene name, biotype from Ensembl
    my $genes = $slice->get_all_Genes();
	
	## the slice can contain multiple genes
	## create 3 arrays to store geneID, name and biotype
    my @geneids;
    my @geneextnames;
    my @genebiotype;

	## loop through $genes and store each geneID, name and biotype
	## in the corresponding array
    while ( my $gene = shift @{$genes} ) {
		push(@geneids, $gene->stable_id());
		push(@geneextnames, $gene->external_name);
		push(@genebiotype, $gene->biotype());
    }	

    ## create arrays to store ENCODE regulatory features:
	## DHS peaks, histone peaks, polymerase peaks and TFBS peaks
    my @dhs;
    my @hist;
    my @pol;
    my @tfbs;

	## create arrays to store the cell types of the ENCODE features
    my @dhs2;
    my @hist2;
    my @pol2;
    my @tfbs2;

	## get regulatory features in this slice
	## loop through all features and store in corresponding array
	my @reg_feats = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};
	foreach my $rf (@reg_feats){
		# print "stable id: ". $rf->stable_id."\n";
		my $rfs = $regfeat_adaptor->fetch_all_by_stable_ID($rf->stable_id);
		foreach my $currentregfeat (@{$rfs}) {

			foreach my $att (@{$currentregfeat->regulatory_attributes()}){ #was rf
				#print_feature($att);
				my $overlap = 0;
				$overlap = ($position>=$att->seq_region_start && $position<=$att->seq_region_end);
				# print $att->display_label, "   snppos: ",$position, " start: ",$att->seq_region_start," end: ",$att->seq_region_end," overlapping with snp: ",$overlap;
				# if ($overlap) {print "   yearh! \n\n"} else {print "\n\n"};
				if ($overlap) {
					my @val = split(/ - /,$att->display_label);
					if    ($val[0] eq 'DNase1')   { push(@dhs,  $val[0])}
					elsif ($val[0] eq 'H3K27ac')  { push(@hist, $val[0])}
					elsif ($val[0] eq 'H3K27me3') { push(@hist, $val[0])}
					elsif ($val[0] eq 'H3K36me3') { push(@hist, $val[0])}
					elsif ($val[0] eq 'H4K91ac')  { push(@hist, $val[0])}
					elsif (substr($val[0],0,3) eq 'H2A')  { push(@hist2, $att->display_label)}
					elsif (substr($val[0],0,3) eq 'H2B')  { push(@hist2, $att->display_label)}
					elsif (substr($val[0],0,3) eq 'H3K')  { push(@hist2, $att->display_label)}
					elsif (substr($val[0],0,3) eq 'H4K')  { push(@hist2, $att->display_label)}
					
					elsif (substr($val[0],0,3) eq 'Pol')  { push(@pol, $val[0])}
					else                          { push(@tfbs, $val[0])};
					if    ($val[0] eq 'DNase1')   { push(@dhs2,  $att->display_label)}
					elsif ($val[0] eq 'H3K27ac')  { push(@hist2, $att->display_label)}
					elsif ($val[0] eq 'H3K27me3') { push(@hist2, $att->display_label)}
					elsif ($val[0] eq 'H3K36me3') { push(@hist2, $att->display_label)}
					elsif (substr($val[0],0,3) eq 'H2A')  { push(@hist2, $att->display_label)}
					elsif (substr($val[0],0,3) eq 'H2B')  { push(@hist2, $att->display_label)}
					elsif (substr($val[0],0,3) eq 'H3K')  { push(@hist2, $att->display_label)}
					elsif (substr($val[0],0,3) eq 'H4K')  { push(@hist2, $att->display_label)}
					
					elsif ($val[0] eq 'H4K91ac')  { push(@hist2, $att->display_label)}
					elsif (substr($val[0],0,3) eq 'Pol')  { push(@pol2, $att->display_label)}
					else                          { push(@tfbs2, $att->display_label)};
				}
			}
		}
	}

	# results are stored in $result
	# if the -t option is set, all results are added to $result separated by tabs
	# if the -t option is not set, all resutls are added to $result separated by semi-colons
	
	# add ENCODE features to $result
	my $result = "";
	if ($opt_t) {
    	$result = join("/",@dhs). "\t".join("/",@hist). "\t".join("/",@pol). "\t".join("/",@tfbs)."\t";
    	$result .=   join("/",@dhs2)."\t".join("/",@hist2)."\t".join("/",@pol2)."\t".join("/",@tfbs2);
	} else {
    	$result = join("/",@dhs). ";".join("/",@hist). ";".join("/",@pol). ";".join("/",@tfbs).";";
    	$result .=   join("/",@dhs2).";".join("/",@hist2).";".join("/",@pol2).";".join("/",@tfbs2);
	}
	
    my @vfs = @{$vfa->fetch_all_by_Slice($slice)};
    my $ref;
    my $reffreq;
    my $alt;
    my $altfreq;

	# add gene annotation to result
	if ($opt_t) {
		$result.="\t".join("/",@geneids)."\t".join("/",@geneextnames)."\t".join("/",@genebiotype);
	} else {
		$result.=";".join("/",@geneids).";".join("/",@geneextnames).";".join("/",@genebiotype);
	}

	# amount of fields in input VCF file
    my $nrfields = @val;
    print OUTFILE $val[0]; # write first field to output file
	# loop over the rest of the fields
	# each time, write a tab and then the field value to the output file
	# only if $i == 7 (the 8th column) we will print $result to the output file
	# first, if populations are detected in the config.txt file,
	# the alleles and frequencies are added to $result
	
    for ($i=1; $i<$nrfields; $i++) { ## a lot of work for only column 7 ;)
	if ($i!=7) {
	    print OUTFILE "\t",$val[$i];
	} else {
	    foreach my $pop (@populations) {
		if (!$opt_t) {
		    $ref     = "${pop}_REF=";
		    $reffreq = "${pop}_REFFREQ=";
		    $alt     = "${pop}_ALT=";
		    $altfreq = "${pop}_ALTFREQ=";
		} else {
		    $ref     = "";
		    $reffreq = "";
		    $alt     = "";
		    $altfreq = "";
		}

		foreach my $vf(@vfs){
		    my @alleles     = @{$vf->variation->get_all_Alleles()};
		    my $id          = $vf->variation_name();
		    my $both        = $vf->allele_string(); 
		    my $refallele   = $vf->ref_allele_string(); 
		    
		    foreach my $allele(@alleles) {
			if($allele->population && $allele->population->name eq $pop ) {
			    # print "currently ref is: ${ref} and the last character is: ", substr($ref, length($ref)-1,1),"\n";
			    # print "currently alt is: ${alt} and the last character is: ", substr($alt, length($alt)-1,1),"\n";
			    
			    if ($allele->allele eq $refallele) {
				if ($ref     eq "" || substr($ref,length($ref)-1,1) eq "=")         { $ref     =            $ref.$allele->allele; }
				else                                                                { $ref     =        $ref."/".$allele->allele; }
				if ($reffreq eq "" || substr($reffreq,length($reffreq)-1,1) eq "=") { $reffreq =     $reffreq.$allele->frequency; }
				else                                                                { $reffreq = $reffreq."/".$allele->frequency; }
			    } else {
				if ($alt     eq "" || substr($alt,length($alt)-1,1) eq "=")         { $alt     =            $alt.$allele->allele; } 
				else                                                                { $alt     =        $alt."/".$allele->allele; }
				if ($altfreq eq "" || substr($altfreq,length($altfreq)-1,1) eq "=") { $altfreq =     $altfreq.$allele->frequency; } 
				else                                                                { $altfreq = $altfreq."/".$allele->frequency; }
			    }
			}
			if($debug == 1 && $allele->population) {
			    print $vf->seq_region_name, "\t", $vf->seq_region_start, "\t", $vf->seq_region_end, "\t",
			    $vf->variation->name, "\t",   $allele->allele, "\t", $vf->allele_string(), " (both) \t", $vf->ref_allele_string(), " (ref) \t",
			    (defined($allele->frequency) ? $allele->frequency : "-"), "\t", $allele->population->name, "\n";
			}
		    }
		}	
		if (!$opt_t) { $result = $result.";".$ref.";".$reffreq.";".$alt.";".$altfreq; }
		if ($opt_t)  { $result = $result."\t".$ref."\t".$reffreq."\t".$alt."\t".$altfreq; }
	    }
		if (!$opt_t) { 
			if ($val[7] eq ".") {
				print OUTFILE "\t",$result; # if input cell was empty, only print $result
			} else {
				print OUTFILE "\t",$val[7], ";", $result;  # if not, print original entry, then semi-colon, then $result
			}
		}	
		if ($opt_t)  { print OUTFILE "\t",$val[7], "\t", $result; }
	    
	}
    }
    print OUTFILE "\n";
}
close(FILE);
close(OUTFILE);

## if($allele->population && $allele->population->name =~ /IIPGA-WEISS-MARTINEZ:D-0/) {
## to test tri-allelic snp   IIPGA-WEISS-MARTINEZ:D-0
