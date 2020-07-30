#!/usr/bin/perl -w
#
use strict;

my $cell_cluster_info_file = $ARGV[0];
my $whitelist = $ARGV[1];
my $bam_file_dir = $ARGV[2];
my $out_file_dir = $ARGV[3];
my $max_distance = $ARGV[4];

printUsage() if(scalar(@ARGV) < 5);

my %barcodes;
my %unplaced_barcodes;

my %cell_info;
my %cluster_info;


system("mkdir -p $out_file_dir");
my @bam_files = glob( $bam_file_dir . '/*bam');


sub main{

  ## Enumerate from white list all possible n distance matches
  open(WHITELIST, "$whitelist") || die("Error reading $whitelist\n");
  while(my $line = <WHITELIST>){
    my @fields = split " ", $line;
    foreach my $s1 (@fields){
      $barcodes{$s1} = 1;
      my @enum = enumerate($s1, $max_distance);
      foreach my $en (@enum){
        my ($string, $poss) = split "\t", $en;
        push(@{$unplaced_barcodes{$string}}, $poss);
      }
    }
  }
  close(WHITELIST);
  
  ## Load the cluster info file
  loadClusterInfo("$cell_cluster_info_file");
  
  ## Open files populate with header lines
  my $header_file = $bam_files[0];
  foreach my $cluster (keys %cluster_info){
    my $filename = $cluster_info{$cluster}->{"file"};
    open($cluster_info{$cluster}->{"filehandle"}, "| samtools view -b - > $filename") || die("cannot open file");
    open(HEADER, "samtools view -H $header_file |") || die("error reading header of bamfile\n");
    while(my $header_line = <HEADER>){
      print { $cluster_info{$cluster}->{"filehandle"} } $header_line;
    }
    close(HEADER);
  }
  
  ## Process one bam at a time 
  foreach my $cur_bam_file (@bam_files){
    splitBamToClusters($cur_bam_file);
  }
  
  ## Close filehandles and sort the bam files, index & make bigWig
  foreach my $cluster (keys %cluster_info){
    close($cluster_info{$cluster}->{"filehandle"});
    my $bam_name = "$out_file_dir/$cluster.bam";
    my $sorted_bam_name = "$out_file_dir/$cluster.sorted.bam";
    my $wig_name = "$out_file_dir/$cluster.bigWig";
    system("samtools sort $bam_name > $sorted_bam_name");
    system("rm $bam_name");
    system("samtools index $sorted_bam_name");
    system("bamCoverage --bam $sorted_bam_name -o $wig_name -bs 50 --normalizeUsing CPM --skipNAs");
  }

}

sub loadClusterInfo{
  my $info_file = shift;
  open(INFILE, "$info_file") || die ("Error reading $info_file\n");
  while(my $line = <INFILE>){
    chomp($line);
    my ($cell, $cluster) = split "\t", $line;
    $cell =~ s/_2/\.2/;
    $cell =~ s/_N/\.N/;
    $cell =~ s/_S/\.S/;
    $cell =~ s/_s/\.s/;
    my ($pool_id, $cell_barcode) = split "_", $cell;
    my $barcode1 = substr($cell_barcode, 0, 8);
    my $barcode2 = substr($cell_barcode, 8, 8);
    my $barcode3 = substr($cell_barcode, 16, 8);

    ## Match the current cell barcodes to whitelist file
    if(!exists($barcodes{$barcode1}) and exists($unplaced_barcodes{$barcode1})){
      my @Poss = @{$unplaced_barcodes{$barcode1}};
      next if(scalar(@Poss) > 1);
      my ($query, $edit) = split ":", $Poss[0];
      $barcode1 = $query;
    }
    if(!exists($barcodes{$barcode2}) and exists($unplaced_barcodes{$barcode2})){
      my @Poss = @{$unplaced_barcodes{$barcode2}};
      next if(scalar(@Poss) > 1);
      my ($query, $edit) = split ":", $Poss[0];
      $barcode2 = $query;
    }
    if(!exists($barcodes{$barcode3}) and exists($unplaced_barcodes{$barcode3})){
      my @Poss = @{$unplaced_barcodes{$barcode3}};
      next if(scalar(@Poss) > 1);
      my ($query, $edit) = split ":", $Poss[0];
      $barcode3 = $query;
    }
    if(exists($barcodes{$barcode1}) and exists($barcodes{$barcode2}) and exists($barcodes{$barcode3})){
      $cell_info{$pool_id}->{$barcode1}->{$barcode2}->{$barcode3} = $cluster;
      $cluster_info{$cluster}->{"file"} = "$out_file_dir/$cluster.bam";
    }
  }
  close(INFILE);
}


sub splitBamToClusters{
  my $bam_file = shift;
  my $m = 0;
  open(BAMFILE, "samtools view $bam_file |") || die("Error reading $bam_file\n");
  while(my $sam_line = <BAMFILE>){
    chomp($sam_line);
    $m++;
    print "$m reads in bam\n" if($m%1000000==0);
    my @fields = split "\t", $sam_line;
    my ($cell, $read_id) = split ":", $fields[0];
    $cell =~ s/_2/\.2/;
    $cell =~ s/_N/\.N/;
    $cell =~ s/_S/\.S/;
    $cell =~ s/_s/\.s/;
    my ($pool_id, $cell_barcode) = split "_", $cell;
    my $barcode1 = substr($cell_barcode, 0, 8);
    my $barcode2 = substr($cell_barcode, 8, 8);
    my $barcode3 = substr($cell_barcode, 16, 8);
    #print $barcode1, ",", $barcode2, ",", $barcode3, "\n";
    if(!exists($barcodes{$barcode1}) and exists($unplaced_barcodes{$barcode1})){
      my @Poss = @{$unplaced_barcodes{$barcode1}};
      next if(scalar(@Poss) > 1);
      my ($query, $edit) = split ":", $Poss[0];
      $barcode1 = $query;
    }
    if(!exists($barcodes{$barcode2}) and exists($unplaced_barcodes{$barcode2})){
      my @Poss = @{$unplaced_barcodes{$barcode2}};
      next if(scalar(@Poss) > 1);
      my ($query, $edit) = split ":", $Poss[0];
      $barcode2 = $query;
    }
    if(!exists($barcodes{$barcode3}) and exists($unplaced_barcodes{$barcode3})){
      my @Poss = @{$unplaced_barcodes{$barcode3}};
      next if(scalar(@Poss) > 1);
      my ($query, $edit) = split ":", $Poss[0];
      $barcode3 = $query;
    }
    if(exists($cell_info{$pool_id})){
      if(exists($cell_info{$pool_id}->{$barcode1})){
        if(exists($cell_info{$pool_id}->{$barcode1}->{$barcode2})){
          if(exists($cell_info{$pool_id}->{$barcode1}->{$barcode2}->{$barcode3})){
            $cell = $pool_id . "_" . $barcode1 . $barcode2 . $barcode3;
            $fields[0] = $cell . ":" . $read_id;
            my $cur_cluster = $cell_info{$pool_id}->{$barcode1}->{$barcode2}->{$barcode3};
            print { $cluster_info{$cur_cluster}->{"filehandle"} } join("\t", @fields), "\n";
          }
        }
      }
    }
  }
  close(BAMFILE);
}


# Return the number of sequence with distance n
#
#
sub enumerate
{
    # $s1 is the sequence
    # $len1 is the length of the sequence
    #
    my ($s1, $n) = @_;
    my ($len1) = length $s1;
    $s1 = uc($s1);
    return 0 if ($len1 == 0 or $n == 0);

    my %mat;
    $mat{$s1} = "NA";
    my @dict = ("A", "T", "C", "G");
    for(my $i = 0; $i < $n; $i++){
      my %enumerated;
      foreach my $string (keys %mat){
        my @edits = split ":", $mat{$string};
        my %mods;
        foreach my $mod (@edits){
          $mods{$mod} = 1;
        }
        if(scalar(@edits) == $i or $i == 0){
          for(my $j = 0; $j < $len1; $j++){
            my $cur_base = substr($s1, $j, 1);
            next if($mods{$j});
            foreach my $alt_base (@dict){
              next if($alt_base eq $cur_base);
              my $cur_seq = substr($string, 0, $j) . $alt_base . substr($string, $j+1, $len1-$j);
              #print $cur_seq, ",", $mat{$string}, ":$j\n";
              if($mat{$string} eq "NA"){
                $mat{$cur_seq} = $j;
              }else{
                $mat{$cur_seq} = $mat{$string}.":".$j;
              }
            }
          }
        }
      }
    }
    my @enum;
    foreach my $string (keys %mat){
      $mat{$string} =~ s/NA://g;
      next if($string eq $s1);
      my @edits = split ":", $mat{$string};
      push(@enum, $string . "\t" . $s1 . ":". scalar(@edits));
    }
    return @enum;
}

sub printUsage{
  print "Usage: makeClusterBigWig.pl <cluster info table> <barcodes whitelist file> <bam file directory> <output directory> <maximum barcode mismatch num> \n";
  exit 1;
}

main;
