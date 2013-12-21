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
    my $seq_region = $feature->slice->seq_region_name();
    my $start      = $feature->start();
    my $end        = $feature->end();
    my $strand     = $feature->strand();

    return sprintf( "%s: %s:%d-%d (%+d)",
		    $stable_id, $seq_region, $start, $end, $strand );
}

my $slice = $slice_adaptor->fetch_by_region( 'chromosome', '1', 100466279 , 100466279 );

my $genes = $slice->get_all_Genes();
while ( my $gene = shift @{$genes} ) {
    my $gstring = feature2string($gene);
    print "$gstring\n";

    my $transcripts = $gene->get_all_Transcripts();
    while ( my $transcript = shift @{$transcripts} ) {
        my $tstring = feature2string($transcript);
        print "\t$tstring biotype: " ,$transcript->biotype(), "\n";

        foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
            my $estring = feature2string($exon);
            print "\t\t$estring\n";
        }
    }
}
