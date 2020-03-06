#!/usr/bin/perl -w 
use strict ;

my $IDfile = shift ;   # ZiYuan Sample ID

my %hashID = () ;
open IN1,"<",$IDfile ;
while (<IN1>){
    chomp ;
    $hashID{"$_"} = 1 ;
}

my %hash = () ;
my %hashALT = () ;
my @files = `ls *.vcf.gz` ;
for my $file (@files){
    chomp($file) ;
    if ($file =~ /(.*?)\.comm\.vcf\.gz/){
        my $id = $1 ;
        my $newName = $id.".comm.vcf" ;
        system("gunzip -c $file >$newName") ;
        if (defined $hashID{"$id"}){
            open IN2,"<",$newName ;
            while (<IN2>){
                chomp ;
                if (! /^#/){
                    my @regions = split ;
                    if ($regions[7] =~ /;DP=(.*?);/){
                        my $dp = $1 ;
                            # select mapping depth >=5 
                        if ($dp >=5){
                            my @values = split(/:/,$regions[9]) ;
                            my $trait = $values[0] ;
                            my $loci = $regions[0].":".$regions[1] ;
                            $hash{"$loci"}->{"REF"} = $regions[3] ;
                            if (! defined $hashALT{"$loci"}->{"ALT"}){
                                $hashALT{"$loci"}->{"ALT"} = $regions[4] ;
                            }
                            $hashALT{"$loci"}->{"$id"} = $trait ;
                        }
                    }
                }
            }   
        }
    }
}

my @sites = keys %hash ;
my @samples = keys %hashID ;
print "Chr","\t","loci","\t","REF" ;
for my $name (@samples){
       print "\t",$name ;
}
print "\n" ;
for my $t (@sites){
    my @tmps = split(/:/,$t) ;
    print $tmps[0],"\t",$tmps[1],"\t",$hash{"$t"}->{"REF"}.$hash{"$t"}->{"REF"} ;
    for my $sample (@samples){
        my $altbase = $hashALT{"$t"}->{"ALT"} ;
        my $altpair = "" ;
        if (defined $hashALT{"$t"}->{"$sample"}){
            my $traitSNP = $hashALT{"$t"}->{"$sample"} ;
            if ($traitSNP eq "0/1"){
                $altpair = $hash{"$t"}->{"REF"}.$altbase ;
            }elsif ($traitSNP eq "1/1"){
                $altpair = $altbase.$altbase ;
            }
            print "\t",$altpair ;
        }else{
            $altpair = $hash{"$t"}->{"REF"}.$hash{"$t"}->{"REF"} ;
            print "\t",$altpair ;
        }
    }
    print "\n" ;
}

