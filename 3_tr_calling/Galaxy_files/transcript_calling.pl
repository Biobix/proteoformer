#!/usr/bin/perl -w

#####################################
##	PROTEOFORMER: deep proteome coverage through ribosome profiling and MS integration
##
##	Copyright (C) 2014 G. Menschaert, J.Crappé, E. Ndah, A. Koch & S. Steyaert
##
##	This program is free software: you can redistribute it and/or modify
##	it under the terms of the GNU General Public License as published by
##	the Free Software Foundation, either version 3 of the License, or
##	(at your option) any later version.
##
##	This program is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##	GNU General Public License for more details.
##
##	You should have received a copy of the GNU General Public License
##	along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
## 	For more (contact) information visit http://www.biobix.be/PROTEOFORMER
#####################################


$|=1;

use strict;
use warnings;
use Getopt::Long;
use v5.10;
use Cwd;

#############
##Command line
#############
#
# perl transcript_calling.pl
#
#

my ($work_dir,$tmp,$in_sqlite,$out_sqlite,$method,$ens_db,$mincount,$no_of_samples,$fdr,$default_score,$cutoff,$alpha,$scripts_folder);

GetOptions(
"work_dir:s"=>\$work_dir,
"tmp:s"=>\$tmp,
"in_sqlite=s"=>\$in_sqlite,
"out_sqlite=s"=>\$out_sqlite,
"method:s"=>\$method,
"ens_db:s"=>\$ens_db,
"min_count:i"=>\$mincount,
"no_of_samples:i"=>\$no_of_samples,
"fdr:f"=>\$fdr,
"default_score:s"=>\$default_score,
"cutoff:f"=>\$cutoff,
"alpha:f"=>\$alpha,
"scripts_folder=s"=>\$scripts_folder
);

my $CWD = getcwd;
if !($work_dir){
    $work_dir = $CWD;
}
if !($tmp){
    $tmp = $work_dir."/tmp/";
}
if !($in_sqlite){
    die "Do not forget to pass the in_sqlite argument!\n";
}
if !($out_sqlite){
    die "Do not forget to pass the out_sqlite argument!\n";
}
if !($method){
    die "Do not forget to pass the method argument!\n";
}
if(lc($method) eq "rule-based"){
    if !($ens_db){
        die "Do not forget to pass the ens_db argument!\n";
    }
} elsif (lc($method) eq "ribozinb"){
    if !($min_count){
        die "Do not forget to pass the min_count argument!\n";
    }
    if !($no_of_samples){
        die "Do not forget to pass the no_of_samples argument!\n";
    }
    if !($fdr){
        die "Do not forget to pass the fdr argument!\n";
    }
    if !($default_score){
        die "Do not forget to pass the default_score argument!\n";
    }
    if !($cutoff){
        die "Do not forget to pass the cutoff argument!\n";
    }
    if !($alpha){
        die "Do not forget to pass the alpha argument!\n";
    }
} else {
    die "Method should be rule-based or ribozinb!\n";
}
if !($scripts_folder){
    die "Do not forget to mention the scripts folder!\n";
}

if (lc($method) eq "rule-based"){
    my $cmd = "perl ".$scripts_folder."/ribo_translation.pl --work_dir ".$work_dir." --tmp ".$tmp." --in_sqlite ".$in_sqlite." --out_sqlite ".$out_sqlite." --ens_db ".$ens_db;
    system($cmd);
} elsif (lc($method) eq "ribozinb"){
    my $cmd = "python ".$scripts_folder."/RiboZINB.py --work_dir ".$work_dir." --tmpfolder ".$tmp." --in_sqlite ".$in_sqlite." --out_sqlite ".$out_sqlite." --mincount ".$mincount." --no_of_samples ".$no_of_samples." --fdr ".$fdr." --default_score ".$default_score." --cutoff ".$cutoff." --alpha ".$alpha;
    system($cmd);
}
