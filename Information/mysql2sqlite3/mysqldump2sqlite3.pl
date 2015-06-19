#!/usr/bin/perl 
#===============================================================================
# ADAPTED FROM:	http://stackoverflow.com/questions/489277/script-to-convert-mysql-dump-sql-file-into-format-that-can-be-imported-into-sqli
# USAGE: 		./mysql2sqlite.pl <MySQL_dumpfile>
# DESCRIPTION: 	Converts MySQL dumpfile to SQLite database
#              	Triggers are not converted
#              	The dump must be done with
#              	> mysqldump --skip-triggers -u [user] --p [database] > dumpfile
# EXAMPLE		> mysqldump --skip-triggers -u you -p ensembl_mus_musculus_core_72_38 gene coord_system exon exon_transcript transcript translation seq_region > ENS_MMU_72.sql				
# REQUIREMENTS: Perl and module SQL::Translator, SQLite
#===============================================================================
use strict;
use warnings;
use Carp;
use English qw( -no_match_vars );
use SQL::Translator;
use 5.012;
use Data::Dumper;

my $file = $ARGV[0];
my $filedb = $file;
$filedb =~ s/\.*[^.]*$/.db/;
if ( -s $filedb ) { 
    say "*** The file already exists < $filedb >. I quit!";
    exit;
}
my @stru;
my @stru_idx;
my @data;
my $tb;

open( my $SQLFILE, "<", $file )
    or croak "Can't open $file: $OS_ERROR";
while (<$SQLFILE>) {
    # not considering lines with comments and lock / unlock / drop
    next if ( /^--/ || /^\/\*/ || /^lock/i || /^unlock/i || /^drop/i );
    # Proces inserts
    if (/^(INSERT.+?)[(]/) {     
        my $ins = $1;            # Capture table name
        s/\\[']/''/g;            # Substitue
        s/[)],[(]/);\n$ins(/g;   # Divide multiple inserts
        push( @data, $_ );
    }
    # Render the structure
    else { push( @stru, $_ ); }
}
close($SQLFILE);

# Ensembl index and table names overlap which isn't supported by SQLite
# Rename indexes to table_name_index	
foreach (@stru) {
 	if($_ =~ /^CREATE TABLE `(.*)`.*/){$tb = $1;}
 	if($_ =~ /(^\s KEY `)(.+?)(`.*)/ ){
 		my $idx = $tb."_".$2;
 		$idx = $1.$idx.$3;	
 		push(@stru_idx, $idx);
 	}elsif(/(^\s UNIQUE KEY `)(.+?)(`.*)/){
 		my $idx = $tb."_".$2;
 		$idx = $1.$idx.$3;
 		push(@stru_idx, $idx);
 	}else{push(@stru_idx,$_);}
 } 

my $strusql = join( '', @stru_idx );
my $datasql = join( '', @data );
#open( my $STRU,   ">", "stru.sql" ); # to verify the results
#open( my $DATA,  ">", "data.sql" );
#print $STRU  $strusql;
#print $DATA  $datasql;

# The conversion, using SQL::Translator
my $translator = SQL::Translator->new(
    no_comments       => 0,
    show_warnings     => 0,
    quote_table_names => 1,
    quote_field_names => 1,
    validate          => 1,
);
my $struout = $translator->translate(
    from => 'MySQL',
    to   => 'SQLite',
    data => \$strusql,
    # filename => $file,
) or croak "Error: " . $translator->error;

# Defines the beginning and end of the INSERTS
my $prgini = "PRAGMA foreign_keys=OFF;\n";
my $traini = "BEGIN TRANSACTION;\n";
my $trafin = "COMMIT;\n";
my $prgfin = "PRAGMA foreign_keys=ON;\n";

# Generates the final SQLITE file
my $sqlout = join( "\n", $struout, $prgini, $traini, $datasql, $trafin, $prgfin);
open( my $FINAL, ">", "/tmp/final.sql" );
print $FINAL $sqlout;

# Create SQLITE database

my $log = "/tmp/sqlite.errlog";
my $command = "sqlite3 $filedb < /tmp/final.sql 2> $log";
system($command) == 0 or die "system $command failed: $?";
if ( -s $log ) { 
    say "*** There was some problem. Check the file < /tmp/sqlite.errlog >"; 
}
else { 
    say "*** Conversion is complete. Check the file < $filedb > "; 
}