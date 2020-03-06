my $file = shift;
open FL,$file;
while(<FL>){
	chomp;
	next if(/^\@/);
	my @tem = split /\t/,$_;
	my $pos = $tem[3];
	my $chain = (!($tem[1] & 0x10)) ? "+":"-";
	my $chr_pos;

	if(!($tem[1] & 0x10)){
#$chr_pos = sprintf "%s:$chain%d",$tem[2],$pos + 17;
		$chr_pos = sprintf "%s:$chain%d",$tem[2],$pos;
	}else{
#$chr_pos = sprintf "%s:$chain%d",$tem[2],$pos + 6;
		$chr_pos = sprintf "%s:$chain%d",$tem[2],$pos;
	}
	next if ($tem[2] eq $chrom && ($pos + 22 <= $end && $pos >= $start));    # in itself
	
	my $str = $tem[12];
	$str =~ s/.*://;
	my @base = split /\d+/,$str;
#my $mm = $#base > 0?$#base:0;
	my @pos = split /[ATCG]/,$str;
	my $site;
	my $foundhit = 1;
	my @mis;
	my @sstr = split //,$tem[9];
	for(1..$#base){
		$site += $pos[$_ -1] + 1;
		$sstr[$site - 1] = $base[$_];
		
		if(!($tem[1] & 0x10)){
			if($site == 22 or $site ==23){
				$foundhit = 0;
			}elsif($site <=20){
				push @mis,$site;
			}
		}else{
			if($site == 2 or $site ==1){
				$foundhit = 0;
			}elsif($site >= 4){
				push @mis, $site;
			}
		}
	}
	if($tem[1] & 0x10){
		for(@sstr){
			tr/ATGC/TACG/;
		}
		@sstr = reverse @sstr;
	}
	my $sstr = join '', @sstr;
	if($foundhit){
		next if($#mis + 1 >=5);
		my $guide_seq = $tem[9];
		my $sstr2 = $sstr;
		$sstr2 =~ s/[^ATCG]//g;
		my $mm = $#mis + 1;
		print "$tem[0]\t$guide_seq\t$sstr\t$mm\t$chr_pos\n";
	}
}

