#modified by Rick Dewey 7/31/12 to produce tab delimited output, and annotate ESP6500, HapMapII+III and all versions of 1000g available in annovar format
#!/usr/bin/perl
use warnings;
use strict;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 504 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2012-05-15 18:05:33 -0700 (Tue, 15 May 2012) $';

our ($verbose, $help, $man);
our ($queryfile, $dbloc);
our ($outfile, $buildver, $step, $checkfile, $remove, $verdbsnp, $ver1000g, $genetype);

GetOptions('verbose|v'=>\$verbose, 'help|h'=>\$help, 'man|m'=>\$man, 'outfile=s'=>\$outfile, 'buildver=s'=>\$buildver, 'step=s'=>\$step, 
	'checkfile!'=>\$checkfile, 'remove'=>\$remove, 'verdbsnp=i'=>\$verdbsnp, 'ver1000g=s'=>\$ver1000g, 'genetype=s'=>\$genetype,) or pod2usage ();
	
$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

my $path = $0;
$path =~ s/[^\\\/]+$//;
$path and $ENV{PATH} = "$path:$ENV{PATH}";		#set up the system executable path to include the path where this program is located in

($queryfile, $dbloc) = @ARGV;

$outfile ||= $queryfile;
$buildver ||= 'hg18';
$buildver eq 'hg18' or $buildver eq 'hg19' or pod2usage ("Error in argument: the --buildver argument can be 'hg18' and 'hg19' only");
not defined $checkfile and $checkfile = 1;

defined $verdbsnp or $verdbsnp = 130;
$genetype ||= 'refgene';
$genetype =~ m/^refgene|knowngene|ensgene|gencodegene$/i or $genetype =~ m/wgEncodeGencode\w+$/ or pod2usage ("Error in argument: the --genetype can be 'refgene', 'knowngene', 'ensgene' or 'gencodegene' only");

my ($file1000g);
if (not defined $ver1000g) {
	if ($buildver eq 'hg18') {
		print STDERR "NOTICE: the --ver1000g argument is set as '1000g' by default\n";
		$ver1000g = '1000g';
	} elsif ($buildver eq 'hg19') {
		print STDERR "NOTICE: the --ver1000g argument is set as '1000g2010nov' by default\n";
		$ver1000g = '1000g2010nov';
	}
}

if ($genetype eq 'gencodegene') {
	if ($buildver eq 'hg18') {
		$genetype = 'wgEncodeGencodeManualV3';
	} elsif ($buildver eq 'hg19') {
		$genetype = 'wgEncodeGencodeManualV4';
	}
}

if ($ver1000g eq '1000g') {
	$file1000g = '2009_04';
	$buildver eq 'hg18' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg18");
} elsif ($ver1000g eq '1000g2010') {
	$file1000g = '2010_03';
	$buildver eq 'hg18' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg18");
} elsif ($ver1000g eq '1000g2010jul') {
	$file1000g = '2010_07';
	$buildver eq 'hg18' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg18");
} elsif ($ver1000g eq '1000g2010nov') {
	$file1000g = '2010_11';
	$buildver eq 'hg19' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg19");
} elsif ($ver1000g eq '1000g2011may') {
	$file1000g = '2011_05';
	$buildver eq 'hg19' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg19");
} elsif ($ver1000g eq '1000g2012feb') {
	$file1000g = '2012_02';
	$buildver eq 'hg19' or pod2usage ("Error in argument: the --ver1000g $ver1000g is supported only for --buildver hg19");
} elsif ($ver1000g =~ m/^1000g(20\d\d)([a-z]{3})_([a-z]+)$/) {
	my %monthhash = ('jan'=>'01', 'feb'=>'02', 'mar'=>'03', 'apr'=>'04', 'may'=>'05', 'jun'=>'06', 'jul'=>'07', 'aug'=>'08', 'sep'=>'09', 'oct'=>'10', 'nov'=>'11', 'dec'=>'12');
	$file1000g = $1 . '_' . $monthhash{$2};
} else {
	pod2usage ("Error in argument: the --ver1000g $ver1000g is not yet supported by this program");
}

my %valistep;
if ($step) {
	my @step = split (/,/, $step);
	for my $i (0 .. @step-1) {
	if ($step[$i] =~ m/^(\d+)-(\d+)$/) {
		for my $nextstep ($1 .. $2) {
		$valistep{$nextstep}++;
		}
	} elsif ($step[$i] =~ m/^(\d+)$/) {
		$valistep{$1}++;
	} else {
		pod2usage ("Error: invalid -step argument ($step) is specified. Please use comma-separated number only (dash line such as 1-5 is accepted)");
	}
	}
} else {
	for my $nextstep (1..35) {
		$valistep{$nextstep}++;
	}
}

$checkfile and checkFileExistence ();


my $sc;

#run step 1
if ($valistep{1}) {
	$sc = "perl $path/annotate_variation.pl -geneanno -buildver $buildver -dbtype $genetype -outfile $outfile -exonsort $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 1 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step2
if ($valistep{2}) {
	if ($buildver eq 'hg18') {
		$sc = "perl $path/annotate_variation.pl -regionanno -dbtype mce44way -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 2 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	} elsif ($buildver eq 'hg19') {
		$sc = "perl $path/annotate_variation.pl -regionanno -dbtype mce46way -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
}

#run step3
if ($valistep{3}) {
	$sc = "perl $path/annotate_variation.pl -regionanno -dbtype segdup -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 3 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}


if ($buildver eq 'hg19') {
	#run step4
	if ($valistep{4}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_ASW.txt -buildver $buildver -outfile $outfile.hg19.hapmap_asw $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 4 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	if ($valistep{5}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_CEU.txt -buildver $buildver -outfile $outfile.hg19.hapmap_ceu $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 5 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	if ($valistep{6}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_CHB.txt -buildver $buildver -outfile $outfile.hg19.hapmap_chb $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 6 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	if ($valistep{7}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_CHD.txt -buildver $buildver -outfile $outfile.hg19.hapmap_chd $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 7 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	if ($valistep{8}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_GIH.txt -buildver $buildver -outfile $outfile.hg19.hapmap_gih $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 8 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	if ($valistep{9}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_JPT.txt -buildver $buildver -outfile $outfile.hg19.hapmap_jpt $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 9 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	if ($valistep{10}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_LWK.txt -buildver $buildver -outfile $outfile.hg19.hapmap_lwk $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 10 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	if ($valistep{11}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_MEX.txt -buildver $buildver -outfile $outfile.hg19.hapmap_mex $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 11 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	if ($valistep{12}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_MKK.txt -buildver $buildver -outfile $outfile.hg19.hapmap_mkk $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 12 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	if ($valistep{13}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_TSI.txt -buildver $buildver -outfile $outfile.hg19.hapmap_tsi $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 13 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	if ($valistep{14}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype generic -genericdbfile hg19.hapmap2and3_YRI.txt -buildver $buildver -outfile $outfile.hg19.hapmap_yri $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 14 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	#run step15
	if ($valistep{15}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype 1000g2010nov_all -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 15 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	#run step16
	if ($valistep{16}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype 1000g2011may_all -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 16 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	#run step17
	if ($valistep{17}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype 1000g2012feb_all -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 17 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	#run step18
	if ($valistep{18}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype 1000g2012apr_all -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 18 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	#run step19
	if ($valistep{19}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype 1000g2012apr_eur -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 19 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
	
	#run step20
	if ($valistep{20}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype 1000g2012apr_afr -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 20 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	#run step21
	if ($valistep{21}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype 1000g2012apr_amr -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 21 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	#run step22
	if ($valistep{22}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype 1000g2012apr_asn -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 22 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	#run step23
	if ($valistep{23}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype cg46 -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 23 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}

	#run step24
	if ($valistep{24}) {
		$sc = "perl $path/annotate_variation.pl -filter -dbtype cg69 -buildver $buildver -outfile $outfile $queryfile $dbloc";
		print STDERR "\nNOTICE: Running step 24 with system command <$sc>\n";
		system ($sc) and die "Error running system command: <$sc>\n";
	}
} 

#run step25
if ($valistep{25}) {
	$sc = "perl $path/annotate_variation.pl -filter -dbtype snp$verdbsnp -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 25 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step26
if ($valistep{26}) {
	$sc = "perl $path/annotate_variation.pl -filter -dbtype avsift -buildver $buildver -sift 0 -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 26 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step27 to step31
if ($valistep{27} or $valistep{28} or $valistep{29} or $valistep{30} or $valistep{31}) {
	$sc = "perl $path/annotate_variation.pl -filter -dbtype ljb_all -buildver $buildver -outfile $outfile $queryfile $dbloc -otherinfo";
	print STDERR "\nNOTICE: Running step 27-31 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step32
if ($valistep{32}) {
	$sc = "perl $path/annotate_variation.pl -filter -dbtype esp6500_all -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 32 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step33
if ($valistep{33}) {
	$sc = "perl $path/annotate_variation.pl -filter -dbtype esp6500_ea -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 33 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

#run step34
if ($valistep{34}) {
	$sc = "perl $path/annotate_variation.pl -filter -dbtype esp6500_aa -buildver $buildver -outfile $outfile $queryfile $dbloc";
	print STDERR "\nNOTICE: Running step 34 with system command <$sc>\n";
	system ($sc) and die "Error running system command: <$sc>\n";
}

open (FUNCTION, "$outfile.variant_function") or die "Error: cannot read from variant function file: $!\n";

if ($valistep{1}) {
	open (STEP1, "$outfile.exonic_variant_function") or die "Error: cannot read from exonic variant function file: $!\n";
}

if ($valistep{2}) {
	if ($buildver eq 'hg18') {
		open (STEP2, "$outfile.hg18_phastConsElements44way") or die "Error: cannot read from mce file: $!\n";
	} 
	elsif ($buildver eq 'hg19') {
		open (STEP2, "$outfile.hg19_phastConsElements46way") or die "Error: cannot read from mce file: $!\n";
	}
}

if ($valistep{3}) {
	open (STEP3, "$outfile.${buildver}_genomicSuperDups") or die "Error: cannot read from segdup file: $!\n";
}

if ($buildver eq 'hg19') {
	if ($valistep{4}) {
		open (STEP4, "$outfile.hg19.hapmap_asw.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_asw.hg19_generic_dropped: $!\n";
	}
	if ($valistep{5}) {
		open (STEP5, "$outfile.hg19.hapmap_ceu.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_ceu.hg19_generic_dropped: $!\n";
	}
	if ($valistep{6}) {
		open (STEP6, "$outfile.hg19.hapmap_chb.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_chb.hg19_generic_dropped: $!\n";
	}
	if ($valistep{7}) {
		open (STEP7, "$outfile.hg19.hapmap_chd.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_chd.hg19_generic_dropped: $!\n";
	}
	if ($valistep{8}) {
		open (STEP8, "$outfile.hg19.hapmap_gih.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_gih.hg19_generic_dropped: $!\n";
	}
	if ($valistep{9}) {
		open (STEP9, "$outfile.hg19.hapmap_jpt.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_jpt.hg19_generic_dropped: $!\n";
	}
	if ($valistep{10}) {
		open (STEP10, "$outfile.hg19.hapmap_lwk.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_lwk.hg19_generic_dropped: $!\n";
	}
	if ($valistep{11}) {
		open (STEP11, "$outfile.hg19.hapmap_mex.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_mex.hg19_generic_dropped: $!\n";
	}
	if ($valistep{12}) {
		open (STEP12, "$outfile.hg19.hapmap_mkk.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_mkk.hg19_generic_dropped: $!\n";
	}
	if ($valistep{13}) {
		open (STEP13, "$outfile.hg19.hapmap_tsi.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_tsi.hg19_generic_dropped: $!\n";
	}
	if ($valistep{14}) {
		open (STEP14, "$outfile.hg19.hapmap_yri.hg19_generic_dropped") or die "Error: cannot read from drop file $outfile.hg19.hapmap_yri.hg19_generic_dropped: $!\n";
	}
	if ($valistep{15}) {
		open (STEP15, "$outfile.hg19_ALL.sites.2010_11_dropped") or die "Error: cannot read from drop file $outfile.hg19_ALL.sites.2010_11_dropped: $!\n";
	}
	if ($valistep{16}) {
		open (STEP16, "$outfile.hg19_ALL.sites.2011_05_dropped") or die "Error: cannot read from drop file $outfile.hg19_ALL.sites.2011_05_dropped: $!\n";
	}
	if ($valistep{17}) {
		open (STEP17, "$outfile.hg19_ALL.sites.2012_02_dropped") or die "Error: cannot read from drop file $outfile.hg19_ALL.sites.2012_02_dropped: $!\n";
	}
	if ($valistep{18}) {
		open (STEP18, "$outfile.hg19_ALL.sites.2012_04_dropped") or die "Error: cannot read from drop file $outfile.hg19_ALL.sites.2012_04_dropped: $!\n";
	}
	if ($valistep{19}) {
		open (STEP19, "$outfile.hg19_EUR.sites.2012_04_dropped") or die "Error: cannot read from drop file $outfile.hg19_EUR.sites.2012_04_dropped: $!\n";
	}
	if ($valistep{20}) {
		open (STEP20, "$outfile.hg19_AFR.sites.2012_04_dropped") or die "Error: cannot read from drop file $outfile.hg19_AFR.sites.2012_04_dropped: $!\n";
	}
	if ($valistep{21}) {
		open (STEP21, "$outfile.hg19_AMR.sites.2012_04_dropped") or die "Error: cannot read from drop file $outfile.hg19_AMR.sites.2012_04_dropped: $!\n";
	}
	if ($valistep{22}) {
		open (STEP22, "$outfile.hg19_ASN.sites.2012_04_dropped") or die "Error: cannot read from drop file $outfile.hg19_ASN.sites.2012_04_dropped: $!\n";
	}
	if ($valistep{23}) {
		open (STEP23, "$outfile.hg19_cg46_dropped") or die "Error: cannot read from drop file $outfile.hg19_cg46_dropped: $!\n";
	}
	if ($valistep{24}) {
		open (STEP24, "$outfile.hg19_cg69_dropped") or die "Error: cannot read from drop file $outfile.hg19_cg69_dropped: $!\n";
	}
} 
else {
	if ($valistep{4}) {
		open (STEP4, "$outfile.hg18_CEU.sites.${file1000g}_dropped") or die "Error: cannot read from drop file $outfile.hg18_CEU.sites.${file1000g}_dropped: $!\n";
	}
	if ($valistep{5}) {
		open (STEP5, "$outfile.hg18_YRI.sites.${file1000g}_dropped") or die "Error: cannot read from drop file $outfile.hg18_YRI.sites.${file1000g}_dropped: $!\n";
	}
	if ($valistep{6}) {
		open (STEP6, "$outfile.hg18_JPTCHB.sites.${file1000g}_dropped") or die "Error: cannot read from drop file $outfile.hg18_JPTCHB.sites.${file1000g}_dropped: $!\n";
	}
}

if ($valistep{25}) {
	open (STEP25, "$outfile.${buildver}_snp${verdbsnp}_dropped") or die "Error: cannot read from snp$verdbsnp drop file: $!\n";
}

if ($valistep{26}) {
	open (STEP26, "$outfile.${buildver}_avsift_dropped") or die "Error: cannot read from avsift drop file: $!\n";
}

if ($valistep{27} or $valistep{28} or $valistep{29} or $valistep{30} or $valistep{31}) {
	open (STEP27, "$outfile.${buildver}_ljb_all_dropped") or die "Error: cannot read from ljb_all drop file: $!\n";
}

if ($valistep{32}) {
	open (STEP32, "$outfile.${buildver}_esp6500_all_dropped") or die "Error: cannot read from esp6500_all drop file: $!\n";
}

if ($valistep{33}) {
	open (STEP33, "$outfile.${buildver}_esp6500_ea_dropped") or die "Error: cannot read from esp6500_ea drop file: $!\n";
}

if ($valistep{34}) {
	open (STEP34, "$outfile.${buildver}_esp6500_aa_dropped") or die "Error: cannot read from esp6500_aa drop file: $!\n";
}

=head1
if ($valistep{11}) {
	open (STEP11, "$outfile.${buildver}_ljb_pp2_dropped") or die "Error: cannot read from ljb_pp2 drop file: $!\n";
}

if ($valistep{12}) {
	open (STEP12, "$outfile.${buildver}_ljb_phylop_dropped") or die "Error: cannot read from ljb_phylop drop file: $!\n";
}

if ($valistep{13}) {
	open (STEP13, "$outfile.${buildver}_ljb_mt_dropped") or die "Error: cannot read from ljb_mt drop file: $!\n";
}

if ($valistep{14}) {
	open (STEP14, "$outfile.${buildver}_ljb_lrt_dropped") or die "Error: cannot read from ljb_lrt drop file: $!\n";
}

if ($valistep{15}) {
	open (STEP15, "$outfile.${buildver}_ljb_gerp++_dropped") or die "Error: cannot read from ljb_gerp++ drop file: $!\n";
}
=cut



my (@allstep);
for my $i (1 .. 3) {
	$allstep[$i] = [];
}

if ($buildver eq 'hg19') {
	if ($valistep{4}) {
		while (<STEP4>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 4 : <$_>\n";
			$allstep[4]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{5}) {
		while (<STEP5>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 5 : <$_>\n";
			$allstep[5]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{6}) {
		while (<STEP6>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 6 : <$_>\n";
			$allstep[6]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{7}) {
		while (<STEP7>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 7 : <$_>\n";
			$allstep[7]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{8}) {
		while (<STEP8>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 8 : <$_>\n";
			$allstep[8]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{9}) {
		while (<STEP9>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 9 : <$_>\n";
			$allstep[9]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{10}) {
		while (<STEP10>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 10 : <$_>\n";
			$allstep[10]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{11}) {
		while (<STEP11>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 11 : <$_>\n";
			$allstep[11]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{12}) {
		while (<STEP12>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 12 : <$_>\n";
			$allstep[12]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{13}) {
		while (<STEP13>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 13 : <$_>\n";
			$allstep[13]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{14}) {
		while (<STEP14>) {
			m/^generic\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 14 : <$_>\n";
			$allstep[14]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{15}) {
		while (<STEP15>) {
			m/^1000g2010nov_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 15 : <$_>\n";
			$allstep[15]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{16}) {
		while (<STEP16>) {
			m/^1000g2011may_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 16 : <$_>\n";
			$allstep[16]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{17}) {
		while (<STEP17>) {
			m/^1000g2012feb_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 17 : <$_>\n";
			$allstep[17]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{18}) {
		while (<STEP18>) {
			m/^1000g2012apr_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 18 : <$_>\n";
			$allstep[18]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{19}) {
		while (<STEP19>) {
			m/^1000g2012apr_eur\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 19 : <$_>\n";
			$allstep[19]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{20}) {
		while (<STEP20>) {
			m/^1000g2012apr_afr\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 20 : <$_>\n";
			$allstep[20]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{21}) {
		while (<STEP21>) {
			m/^1000g2012apr_amr\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 21 : <$_>\n";
			$allstep[21]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{22}) {
		while (<STEP22>) {
			m/^1000g2012apr_asn\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 22 : <$_>\n";
			$allstep[22]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{23}) {
		while (<STEP23>) {
			m/^cg46\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 23 : <$_>\n";
			$allstep[23]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	if ($valistep{24}) {
		while (<STEP24>) {
			m/^cg69\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 24 : <$_>\n";
			$allstep[24]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
} 
else {
	if ($valistep{4}) {
		while (<STEP4>) {
			m/^1000g\w*_ceu\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 4 : <$_>\n";
			$allstep[4]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}

	if ($valistep{5}) {
		while (<STEP5>) {
			m/^1000g\w*_yri\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 5 : <$_>\n";
			$allstep[5]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
	
	if ($valistep{6}) {
		while (<STEP6>) {
			m/^1000g\w*_jptchb\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 6 : <$_>\n";
			$allstep[6]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
		}
	}
}

if ($valistep{25}) {
	while (<STEP25>) {
		m/^snp\d+\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 25 : <$_>\n";
		$allstep[25]->{$2} = $1;
	}
}

if ($valistep{26}) {
	while (<STEP26>) {
		m/^avsift\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 26 : <$_>\n";
		$allstep[26]->{$2} = $1;
	}
}

if ($valistep{27} or $valistep{28} or $valistep{29} or $valistep{30} or $valistep{31}) {
	while (<STEP27>) {
		m/^ljb_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 27 : <$_>\n";
		$allstep[27]->{$2} = $1;
	}
}

if ($valistep{32}) {
	while (<STEP32>) {
		m/^esp6500_all\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 32 : <$_>\n";
		$allstep[32]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
	}
}

if ($valistep{33}) {
	while (<STEP33>) {
		m/^esp6500_ea\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 33 : <$_>\n";
		$allstep[33]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
	}
}

if ($valistep{34}) {
	while (<STEP34>) {
		m/^esp6500_aa\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 34 : <$_>\n";
		$allstep[34]->{$2} = ($1>0.01)?sprintf("%.2f", $1):$1;
	}
}

=head1
if ($valistep{11}) {
	while (<STEP11>) {
		m/^ljb_pp2\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 11 : <$_>\n";
		$allstep[11]->{$2} = $1;
	}
}

if ($valistep{12}) {
	while (<STEP12>) {
		m/^ljb_phylop\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 12 : <$_>\n";
		$allstep[12]->{$2} = $1;
	}
}

if ($valistep{13}) {
	while (<STEP13>) {
		m/^ljb_mt\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 13 : <$_>\n";
		$allstep[13]->{$2} = $1;
	}
}

if ($valistep{14}) {
	while (<STEP14>) {
		m/^ljb_lrt\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 14 : <$_>\n";
		$allstep[14]->{$2} = $1;
	}
}

if ($valistep{15}) {
	while (<STEP15>) {
		m/^ljb_gerp\+\+\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 15 : <$_>\n";
		$allstep[15]->{$2} = $1;
	}
}
=cut

print STDERR "NOTICE: Finished loading filterstep database file\n";

open (OUT, ">$outfile.genome_summary.txt") or die "Error: cannot write to output file: $!\n";
open (OUTE, ">$outfile.exome_summary.txt") or die "Error: cannot write to output file: $!\n";

if ($buildver eq 'hg19') {
	print OUT join ("\t", qw/Func Gene ExonicFunc AAChange Conserved SegDup/, "ESP6500_ALL", "ESP6500_EA", "ESP6500_AA", "hapmap2and3_ASW", "hapmap2and3_CEU", "hapmap2and3_CHB", "hapmap2and3_CHD", "hapmap2and3_GIH", "hapmap2and3_JPT", "hapmap2and3_LWK", "hapmap2and3_MEX", "hapmap2and3_MKK", "hapmap2and3_TSI", "hapmap2and3_YRI", "1000g2010nov_ALL", "1000g2011may_ALL", "1000g2012feb_ALL", "1000g2012apr_ALL", "1000g2012apr_EUR", "1000g2012apr_AFR", "1000g2012apr_AMR", "1000g2012apr_ASN", "cg46", "cg69", "dbSNP$verdbsnp", qw/AVSIFT LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++ Chr Start End Ref Obs Otherinfo/), "\n";
	print OUTE join ("\t", qw/Func Gene ExonicFunc AAChange Conserved SegDup/, "ESP6500_ALL", "ESP6500_EA", "ESP6500_AA", "hapmap2and3_ASW", "hapmap2and3_CEU", "hapmap2and3_CHB", "hapmap2and3_CHD", "hapmap2and3_GIH", "hapmap2and3_JPT", "hapmap2and3_LWK", "hapmap2and3_MEX", "hapmap2and3_MKK", "hapmap2and3_TSI", "hapmap2and3_YRI", "1000g2010nov_ALL", "1000g2011may_ALL", "1000g2012feb_ALL", "1000g2012apr_ALL", "1000g2012apr_EUR", "1000g2012apr_AFR", "1000g2012apr_AMR", "1000g2012apr_ASN", "cg46", "cg69", "dbSNP$verdbsnp", qw/AVSIFT LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++ Chr Start End Ref Obs Otherinfo/), "\n";
} else {
	print OUT join ("\t", qw/Func Gene ExonicFunc AAChange Conserved SegDup/, "ESP6500_ALL", "${ver1000g}_CEU,${ver1000g}_YRI,${ver1000g}_JPTCHB", "dbSNP$verdbsnp", qw/AVSIFT LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++ Chr Start End Ref Obs Otherinfo/), "\n";
	print OUTE join ("\t", qw/Func Gene ExonicFunc AAChange Conserved SegDup/, "ESP6500_ALL", "${ver1000g}_CEU,${ver1000g}_YRI,${ver1000g}_JPTCHB", "dbSNP$verdbsnp", qw/AVSIFT LJB_PhyloP LJB_PhyloP_Pred LJB_SIFT LJB_SIFT_Pred LJB_PolyPhen2 LJB_PolyPhen2_Pred LJB_LRT LJB_LRT_Pred LJB_MutationTaster LJB_MutationTaster_Pred LJB_GERP++ Chr Start End Ref Obs Otherinfo/), "\n";
}

while (<FUNCTION>) {
	s/[\r\n]+$//;
	m/^(\S+)\t([^\t]+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)(.*)/ or die "Error: invalid record found in annovar outputfile: <$_>\n";
	my ($function, $gene, $varstring, $otherinfo) = ($1, $2, $3, $4||'');
	my $exonic;
	if ($function =~ m/\bsplicing\b/ or $function =~ m/\bexonic\b/) {
		$exonic = 1;
	}
	print OUT "$function"."\t"."$gene";
	$exonic and print OUTE "$function"."\t"."$gene";
	
	if (not @{$allstep[1]}) {
		if ($valistep{1}) {
			if (defined ($_ = <STEP1>)) {
				m/^line\d+\t([^\t]+)\t(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 2: <$_>\n";
				my ($efun, $aachange, $varstring) = ($1, $2, $3);
				my @aachange = split (/:|,/, $aachange);
				if (@aachange >= 5) {
					push @{$allstep[1]}, $varstring, $efun, "$aachange[1]:$aachange[3]:$aachange[4]";
				} else {
					push @{$allstep[1]}, $varstring, $efun, $aachange;		#aachange could be "UNKNOWN"
				}
			}
		}
	}
	
	if (not @{$allstep[2]}) {
		if ($valistep{2}) {
			if (defined ($_ = <STEP2>)) {
				m/^mce\d+way\tScore=(\S+)\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 2: <$_>\n";
				push @{$allstep[2]}, $2, $1;
			}
		}
	}
	
	if (not @{$allstep[3]}) {
		if ($valistep{3}) {
			if (defined ($_ = <STEP3>)) {
				m/^segdup\tScore=(\S+);\S+\t(\S+\s+\S+\s+\S+\s+\S+\s+\S+)/ or die "Error: invalid record found in annovar outputfile 3 : <$_>\n";
				push @{$allstep[3]}, $2, ($1>0.01)?sprintf("%.2f", $1):$1;
			}
		}
	}
	
	for my $i (1 .. 3) {			
		my $curstep = $allstep[$i];
		if (@$curstep and $curstep->[0] eq $varstring) {
			if ($i == 1) {
				print OUT "\t"."$curstep->[1]"."\t"."$curstep->[2]";
				$exonic and print OUTE "\t"."$curstep->[1]"."\t"."$curstep->[2]";
			} else {
				print OUT "\t"."$curstep->[1]";
				$exonic and print OUTE "\t"."$curstep->[1]";
			}
			@$curstep = ();
		}
		else {
			if ($i == 1) {
				print OUT "\t"."\t";
				$exonic and print OUTE "\t"."\t";
			} else {
				print OUT "\t";
				$exonic and print OUTE "\t";
			}
		}
	}
	
	if (defined $allstep[32]->{$varstring}) {
		print OUT "\t"."$allstep[32]->{$varstring}";
		$exonic and print OUTE "\t"."$allstep[32]->{$varstring}";
	} else {
		print OUT "\t";
		$exonic and print OUTE "\t";
	}

	if (defined $allstep[33]->{$varstring}) {
		print OUT "\t"."$allstep[33]->{$varstring}";
		$exonic and print OUTE "\t"."$allstep[33]->{$varstring}";
	} else {
		print OUT "\t";
		$exonic and print OUTE "\t";
	}

	if (defined $allstep[34]->{$varstring}) {
		print OUT "\t"."$allstep[34]->{$varstring}";
		$exonic and print OUTE "\t"."$allstep[34]->{$varstring}";
	} else {
		print OUT "\t";
		$exonic and print OUTE "\t";
	}
	
	for my $i (4 .. 26) {		
		if (defined $allstep[$i]->{$varstring}) {
			print OUT "\t"."$allstep[$i]->{$varstring}";
			$exonic and print OUTE "\t"."$allstep[$i]->{$varstring}";
		} else {
			print OUT "\t";
			$exonic and print OUTE "\t";
		}
	}

	#step27 already includes step 27 through step 31
	if (defined $allstep[27]->{$varstring}) {
		print OUT "\t", join("\t", split(/,/, $allstep[27]->{$varstring}));
		$exonic and print OUTE "\t", join("\t", split(/,/, $allstep[27]->{$varstring}));
	} else {
		print OUT "\t"."\t"."\t"."\t"."\t"."\t"."\t"."\t"."\t"."\t"."\t";
		$exonic and print OUTE "\t"."\t"."\t"."\t"."\t"."\t"."\t"."\t"."\t"."\t"."\t";
	}
		
	my @varstring = split (/\s+/, $varstring);
	$otherinfo =~ s/^\s+//;
	my @otherinfo = split (/\t/, $otherinfo);
	for my $i (0 .. @otherinfo-1) {
		$otherinfo[$i] = "$otherinfo[$i]";
	}

	print OUT "\t", join ("\t", @varstring), "\t", join ("\t", @otherinfo), "\n";
	$exonic and print OUTE "\t", join ("\t", @varstring), "\t", join ("\t", @otherinfo), "\n";
}

print STDERR "NOTICE: Final whole-genome summary was written to $outfile.genome_summary.txt file\n";
print STDERR "NOTICE: Final whole-exome summary was written to $outfile.exome_summary.txt file\n";

if ($remove) {
	unlink ("$outfile.variant_function", "$outfile.exonic_variant_function", "$outfile.hg18_phastConsElements44way", "$outfile.hg19_phastConsElements46way", 
		"$outfile.hg19_genomicSuperDups", 
		"$outfile.hg19_esp6500_all_dropped",
		"$outfile.hg19_esp6500_ea_dropped",
		"$outfile.hg19_esp6500_aa_dropped",
		"$outfile.hg19_hapmap2and3_asw.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_ceu.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_chb.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_chd.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_gih.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_jpt.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_lwk.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_mex.hg19_generic_dropped",	
		"$outfile.hg19_hapmap2and3_mkk.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_tsi.hg19_generic_dropped",
		"$outfile.hg19_hapmap2and3_yri.hg19_generic_dropped",
		"$outfile.hg19_ALL.sites.2010_11_dropped",
		"$outfile.hg19_ALL.sites.2011_05_dropped",
		"$outfile.hg19_ALL.sites.2012_02_dropped",
		"$outfile.hg19_ALL.sites.2012_04_dropped",
		"$outfile.hg19_EUR.sites.2012_04_dropped",
		"$outfile.hg19_AFR.sites.2012_04_dropped",
		"$outfile.hg19_AMR.sites.2012_04_dropped",
		"$outfile.hg19_ASN.sites.2012_04_dropped",
		"$outfile.hg19_esp6500_all_dropped",
		"$outfile.hg19_esp6500_ea_dropped",
		"$outfile.hg19_esp6500_aa_dropped",
		"$outfile.hg19_hapmap2and3_asw.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_ceu.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_chb.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_chd.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_gih.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_jpt.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_lwk.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_mex.hg19_generic_filtered",	
		"$outfile.hg19_hapmap2and3_mkk.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_tsi.hg19_generic_filtered",
		"$outfile.hg19_hapmap2and3_yri.hg19_generic_filtered",
		"$outfile.hg19_ALL.sites.2010_11_filtered",
		"$outfile.hg19_ALL.sites.2011_05_filtered",
		"$outfile.hg19_ALL.sites.2012_02_filtered",
		"$outfile.hg19_ALL.sites.2012_04_filtered",
		"$outfile.hg19_EUR.sites.2012_04_filtered",
		"$outfile.hg19_AFR.sites.2012_04_filtered",
		"$outfile.hg19_AMR.sites.2012_04_filtered",
		"$outfile.hg19_ASN.sites.2012_04_filtered", 
		"$outfile.hg19_esp6500_all_filtered",
		"$outfile.hg19_esp6500_ea_filtered",
		"$outfile.hg19_esp6500_aa_filtered",
		"$outfile.hg18_CEU.sites.${file1000g}_dropped", "$outfile.hg18_YRI.sites.${file1000g}_dropped", "$outfile.hg18_JPTCHB.sites.${file1000g}_dropped", 
		"$outfile.${buildver}_snp${verdbsnp}_dropped", "$outfile.${buildver}_avsift_dropped", "$outfile.${buildver}_ljb_all_dropped");
}

sub checkFileExistence {
	my @file = ("${buildver}_refGene.txt", "${buildver}_refLink.txt", "${buildver}_refGeneMrna.fa", "${buildver}_genomicSuperDups.txt", 
		"${buildver}_snp$verdbsnp.txt", "${buildver}_avsift.txt", "${buildver}_ljb_all.txt",  "${buildver}_esp6500_all.txt");
	if ($buildver eq 'hg18') {
		push @file, "${buildver}_phastConsElements44way.txt";
		push @file, "${buildver}_CEU.sites.${file1000g}.txt", "${buildver}_YRI.sites.${file1000g}.txt", "${buildver}_JPTCHB.sites.${file1000g}.txt";
	} elsif ($buildver eq 'hg19') {
		push @file, "${buildver}_phastConsElements46way.txt";
		push @file, "${buildver}_ALL.sites.2012_04.txt", "${buildver}_EUR.sites.2012_04.txt", "${buildver}_AFR.sites.2012_04.txt", "${buildver}_AMR.sites.2012_04.txt", "${buildver}_ASN.sites.2012_04.txt";
	}
	for my $i (0 .. @file-1) {
		my $dbfile = File::Spec->catfile ($dbloc, $file[$i]);
		-f $dbfile or die "Error: the required database file $dbfile does not exist. Please download it via -downdb argument by annotate_variation.pl.\n";
	}
}

=head1 SYNOPSIS

 summarize_annovar.pl [arguments] <query-file> <database-location>

 Optional arguments:
		-h, --help			print help message
		-m, --man			print complete documentation
		-v, --verbose			use verbose output
		    --outfile <string>		output file name prefix
		    --buildver <string>		genome build version (default: hg18)
		    --remove			remove all temporary files
		    --verdbsnp <int>		dbSNP version to use (default: 130)
		    --ver1000g <string>		1000G version (default: 1000g2010nov)
		    --genetype <string>		gene definition can be refgene (default), knowngene, ensgene
		    --checkfile			check existence of database file (default: ON)

 Function: automatically run a pipeline on a list of variants and summarize 
 their functional effects in a comma-delimited file, to be opened by Excel for 
 manual filtering
 
 Example: summarize_annovar.pl ex2.human humandb/
 
 Version: $LastChangedDate: 2012-05-15 18:05:33 -0700 (Tue, 15 May 2012) $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--outfile>

the prefix of output file names

=item B<--buildver>

specify the genome build version

=item B<--remove>

remove all temporary files. By default, all temporary files will be kept for 
user inspection, but this will easily clutter the directory.

=item B<--verdbsnp>

version of dbSNP to use in annotation. By default, 130 is used.

=item B<--ver1000g>

version of 1000 Genomes Project dataset to use in annotation. By default, 1000g2010nov is used.

=item B<--genetype>

gene definition systems to use, including refgene (default), knowngene, ensgene

=item B<--checkfile>

check to make sure that database files exist, before executing the current 
program.

=back

=head1 DESCRIPTION

ANNOVAR is a software tool that can be used to functionally annotate a list of 
genetic variants, possibly generated from next-generation sequencing 
experiments. For example, given a whole-genome resequencing data set for a human 
with specific diseases, typically around 3 million SNPs and around half million 
insertions/deletions will be identified. Given this massive amounts of data (and 
candidate disease- causing variants), it is necessary to have a fast algorithm 
that scans the data and identify a prioritized subset of variants that are most 
likely functional for follow-up Sanger sequencing studies and functional assays.

summarize_annovar is a script that automate some routines in ANNOVAR and 
generates an Excel-compatible file for users to manually browse and filter.

ANNOVAR is freely available to the community for non-commercial use. For 
questions or comments, please contact kai@openbioinformatics.org.


=cut
