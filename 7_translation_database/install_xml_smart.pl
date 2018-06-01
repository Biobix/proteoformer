$|=1;
#!/usr/bin/perl

use strict;
use warnings;

my $install_dir = $INC[0];
print $install_dir;

system("wget http://search.cpan.org/CPAN/authors/id/G/GM/GMPASSOS/Object-MultiType-0.05.tar.gz");
system("gunzip Object-MultiType-0.05.tar.gz");
system("tar -xvf Object-MultiType-0.05.tar");
system("rm -rf Object-MultiType-0.05.tar");
system("mv Object-MultiType-0.05 ".$install_dir);

system("wget http://search.cpan.org/CPAN/authors/id/T/TM/TMHARISH/XML-Smart-1.79.tar.gz");
system("gunzip XML-Smart-1.79.tar.gz");
system("tar -xvf XML-Smart-1.79.tar");
system("rm -rf XML-Smart-1.79.tar");
system("mv XML-Smart-1.79 ".$install_dir);

chdir $install_dir."/Object-MultiType-0.05";
system("perl Makefile.PL");
system("make");
system("make test");
system("make install");

chdir $install_dir."/XML-Smart-1.79";
system("perl Makefile.PL");
system("make");
system("make test");
system("make install");
