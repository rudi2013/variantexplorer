use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
    #The Registry automatically picks port 5306 for ensembl db
    #-verbose => 1, #specificy verbose to see exactly what it loaded
    );

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );
my $tr_adaptor    = $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
my $daf_adaptor   = $registry->get_adaptor( 'Human', 'Core', 'DnaAlignFeature' );

sub feature2string
{
    my $feature = shift;

    my $stable_id  = $feature->stable_id();
    my $extname    = $feature->external_name;
    my $type       = $feature->biotype();
    my $desc       = defined($feature->description()) ? $feature->description() : "";
    my $seq_region = $feature->slice->seq_region_name();
    my $start      = $feature->start();
    my $end        = $feature->end();
    my $strand     = $feature->strand();

    return sprintf( "%s\t%s\t%s\t%s",
		    $stable_id, $extname, $type, $desc );
}

## NOD2 ENSG00000167207
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', '16', 50727514, 50766988 );

my $genes = $slice->get_all_Genes();
while ( my $gene = shift @{$genes} ) {
    my $gstring = feature2string($gene);
    print "$gstring\n";
}

## lincRNA
$slice = $slice_adaptor->fetch_by_region( 'chromosome', '1', 1365990 , 1365990 );

$genes = $slice->get_all_Genes();
while ( my $gene = shift @{$genes} ) {
    my $gstring = feature2string($gene);
    print "$gstring\n";
}
