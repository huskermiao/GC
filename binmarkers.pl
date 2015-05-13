#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

#--Preset--#
my $author = "Leiting Li";
my $version = "20150513";

sub message {local $\ = "\n"; print STDERR @_}
sub usage{
    print <<EOF;

perl binmarkers.pl [OPTIONS] <MARKER_MATRIX>

  MARKER_MATRIX  Scafold/contig and position information
                 as marker name, '-' or '_' as seperator,

  OPTIONS
  
  -t,--threshold  
        Maximum difference allowed within a block,
        default: 5

  -e,--heterozygous 
        Specify heterogygous genotypes
        Default: lm,np,hk,h
        
  -o,--homologous
        Specify homologous genotypes
        Default: ll,nn,hh,kk,a,b

  -m,--missing
        Specify missing genotypes
        Default: --,..,-,.,u
        
  -h, --help      
        print this usage message
  
EOF
    exit;
}

my $infile;
my $threshold = 5;
my $help;
my @heterozygous_genotype = qw/lm np hk h/;
my @homologous_genotype = qw/ll nn hh kk a b/;
my @missing_genotype = qw/-- .. - . u/;

GetOptions(
  "t|threshold=s" => \$threshold,
  "h|help" => \$help,
  "e|heterozygous=s" => \@heterozygous_genotype,
  "o|homologous=s" => \@homologous_genotype,
  "m|missing=s" => \@missing_genotype
) or die("Error in command line arguments");
@heterozygous_genotype = split(/,/,join(',',@heterozygous_genotype));
@homologous_genotype = split(/,/,join(',',@homologous_genotype));
@missing_genotype = split(/,/,join(',',@missing_genotype));

usage if $help;

$infile = shift @ARGV;
usage unless $infile;

my @valid_genotype = (@homologous_genotype, @heterozygous_genotype);
my @all_genotypes = (@valid_genotype, @missing_genotype);
my %is_homologous = map{$_ => 1}@homologous_genotype;
my %is_heterozygous = map{$_ => 1}@heterozygous_genotype;
my %is_valid = map{$_ => 1}@valid_genotype;
my %is_missing = map{$_ => 1}@missing_genotype;
my %is_genotype = map{$_ => 1}@all_genotypes;

my %matrix;
my @marker_index;
my $num_of_markers = 0;
my %bin_marker;
my $num_of_rand_select = 0;
my $num_of_filled_missing = 0;
my $num_of_heterozygous_select = 0;
my $num_of_countif_select = 0;

my $log_file = "$infile.bin_markers_$threshold.log";
my $bin_markers_file = "$infile.bin_markers_$threshold.map";
my $manual_checking_file = "$infile.bin_markers_$threshold.check.txt";

my $status_memory = 0;
sub print_status{
    my ($current, $total) = @_;
    if($current / $total * 10 > $status_memory){
        $status_memory++;
        message "${status_memory}0%";
    }
}

my %bin_marker_threshold;
my $total_diff_cal_times;
sub difference{
    my ($block_id, $a, $b) = @_;
    $total_diff_cal_times++;
    my @a = @{$matrix{$a}->{array}};
    my @b = @{$matrix{$b}->{array}};
    for my $index (1..$#a){
         my $element_a = $a[$index];
         my $element_b = $b[$index];
         if($is_missing{$element_a} or 
             $is_missing{$element_b} or 
             $element_a ne $element_b){
            $bin_marker_threshold{$block_id}->{$index}++;
         }
    }
    my $difference = keys %{$bin_marker_threshold{$block_id}};
    #message "Block $block_id: Difference for $a and $b is $difference";
    return $difference;
}

sub with_title{
    my $array = shift;
    if($is_genotype{$array->[1]}){
        return 0
    }else{
        return 1
    }
}

sub checking_genotypes{
    my $array = shift;
    my @genotypes = @{$array}[1..$#{$array}];
    for (@genotypes){
        die "CAUTION: $_ is not an valid genotype!" 
            unless $is_genotype{$_};
    }
}

sub count_missing{
    my $array = shift;
    my @genotypes = @{$array}[1..$#{$array}];
    my $missing = 0;
    for(@genotypes){
        $missing++ if $is_missing{$_}
    }
    return $missing;
}

open my $log_fh, "> $log_file" or die;
sub print_log{
    print $log_fh @_,"\n";
}

my $title;
sub load_marker_matrix{
    my $file = shift;
    open my $fh, $file or die;
    while(my $marker = <$fh>){
        chomp $marker;
        my $ref = [split /\s+/, $marker];
        my ($scaffold, $position) = split /[\-_]/, $ref->[0];
        if($num_of_markers == 0 and with_title $ref){
            $title = "$marker\n";
            next;
        }
        checking_genotypes $ref;
        my $index = ++$num_of_markers;
        push @marker_index, $index;
        $matrix{$index}->{array} = $ref;
        $matrix{$index}->{missing} = count_missing $ref;
        $matrix{$index}->{scaffold} = $scaffold;
        $matrix{$index}->{position} = $position;
    }
}

sub random_select{
    my $n = scalar(@_);
    die unless $n > 0;
    my $random_index = int(rand($n));
    return $_[$random_index];
}

sub equal_case_select{
    for(@_){
        if($is_heterozygous{$_}){
            $num_of_heterozygous_select++;
            return $_;
        }
    }
    $num_of_rand_select++;
    return random_select(@_);
}

sub select_genotype{
    my @valid_genotypes = @{$_[0]};
    my %countif = %{$_[1]};
    my ($first, $second, $third) = 
        sort{$countif{$b} <=> $countif{$a}} @valid_genotypes;
    
    if(@valid_genotypes == 2){
        if($countif{$first} > $countif{$second}){
            $num_of_countif_select++;
            return $first;
        }else{
            return equal_case_select($first, $second);
        }
    }elsif(@valid_genotypes == 3){
        if($countif{$first} > $countif{$second}){
            $num_of_countif_select++;
            return $first;
        }elsif($countif{$second} > $countif{$third}){
            return equal_case_select($first, $second);
        }else{
            return equal_case_select($first, $second, $third);
        } 
    }else{die @valid_genotypes}
}

sub judge_genotype{
    my @genotypes = @{$_[0]};
    my %countif;
    map{$countif{$_}++}@genotypes;
    my @all_genotypes = keys %countif;
    my @valid_genotypes = grep{$is_valid{$_}}@all_genotypes;
    my @missing_genotypes = grep{$is_missing{$_}}@all_genotypes;
    my $result;
    if(@all_genotypes == 1){
        return $all_genotypes[0];
    }elsif(@valid_genotypes == 1){
        $num_of_filled_missing++ if @missing_genotypes;
        return $valid_genotypes[0];
    }elsif(@valid_genotypes == 2 or @valid_genotypes == 3){
        $num_of_filled_missing++ if @missing_genotypes;
        return select_genotype(\@valid_genotypes, \%countif);
    }else{die "More than three genotypes? @valid_genotypes"}
}

sub consensus_marker{
    my $block_id = shift;
    my @marker_indexes = @_;
    my @consensus_marker;

    my $max_array_index = $#{$matrix{$marker_indexes[0]}->{array}};
    for(my $i = 1; $i <= $max_array_index; $i++){
        my @genotypes = map{$matrix{$_}->{array}->[$i]}@marker_indexes;
        push @consensus_marker, judge_genotype(\@genotypes);
    }
    return @consensus_marker;
}

my %block_contents;
sub cluster_markers{
    my @sorted_marker_index = sort
        {$matrix{$a}->{scaffold} cmp $matrix{$b}->{scaffold} or
         $matrix{$a}->{position} <=> $matrix{$b}->{position}
        }keys %matrix;
    my $block_id = 0;
    for (my $i = 0; $i <= $#sorted_marker_index; $i++){
        my $marker_index = $sorted_marker_index[$i];
        print_status($i, $num_of_markers);
        my $start_scaffold = $matrix{$marker_index}->{scaffold};
        my $end_scaffold;
        my $start_position = $matrix{$marker_index}->{position};
        my $end_position = $start_position;

        #-- Start dynanmic searching 
        $block_id++;
        my @block;
        push @block, $marker_index;
        LABEL: for(my $j = $i + 1; $j <= $#sorted_marker_index; $j++){
            my $new_marker_index = $sorted_marker_index[$j];
            $end_scaffold = $matrix{$new_marker_index}->{scaffold};
            
            # Stop searching with these conditions
            last LABEL if $end_scaffold ne $start_scaffold;
            for my $index_in_block (@block){
                my $d = difference($block_id, $index_in_block, $new_marker_index);
                last LABEL if $d > $threshold;
            }
            
            # Add a marker to a block
            push @block, $new_marker_index;
            $i++;
            $end_position = $matrix{$new_marker_index}->{position};
        }
        #-- End dynamic searching
        
        map{$matrix{$_}->{cluster} = $block_id}@block;
        $block_contents{$block_id} = [@block];
        my $block_size = scalar(@block);
        my $pos_info = $block_size == 1 ? 
                           $start_position : 
                           "$start_position-$end_position";
        my $marker_name = "bin_${start_scaffold}_$pos_info($block_size)";
        $bin_marker{$block_id} = [$marker_name, consensus_marker($block_id, @block)];

        print_log("$marker_name: ", join(", ", map{$matrix{$_}->{array}->[0]}@block));
    }
    return $block_id;
}

# use $out_file and %bin_marker
sub print_bin_markers{
    my $out_file = shift;
    open my $out_fh, ">", $out_file or die;
    print $out_fh $title if $title;
    for(sort{$a <=> $b}(keys %bin_marker)){
        print $out_fh join("\t", @{$bin_marker{$_}}), "\n";
    }
    close $out_fh;
}

# use %matrix and 
sub print_marker_matrix{
    my $out_file = shift;
    open my $out_fh, ">", $out_file or die;
    print $out_fh $title if $title;
    for my $block_id (sort {$a <=> $b} keys %block_contents){
        for my $index (@{$block_contents{$block_id}}){
            my $bin_marker_name = $bin_marker{$block_id}->[0];
            print $out_fh "$bin_marker_name|",join("\t", @{$matrix{$index}->{array}}), "\n";
        }
        print $out_fh join("\t", @{$bin_marker{$block_id}}), "\n";
    }
    close $out_fh;
}

sub hr{message '-' x 60}

sub main{
    message "Loading marker matrix ...";
    load_marker_matrix($infile);
    message "Done! $num_of_markers markers";

    message "Start clustering ...";
    my $num_of_bin_markers = cluster_markers;
    message "Done! $num_of_bin_markers clusters!";

    message "Print original markers append bin marker name";
    print_marker_matrix($manual_checking_file);
    message "Done!";

    message "Print bin markers ...";
    print_bin_markers($bin_markers_file);
    message "Done!";

    hr;
    message "Some log information is in file: $log_file";
    message "Marker matrix for manual checking is in file: $manual_checking_file";
    message "Final bin markers for further analysis if in file: $bin_markers_file";
    hr;
    message "Total number of markers: $num_of_markers";
    message "Bin markers: $num_of_bin_markers";
    message "$num_of_countif_select genotypes were determined based on which present the most";
    message "$num_of_rand_select genotypes were random determined!";
    message "$num_of_filled_missing missing genotypes were filled";
    message "Determination based on heterozygous markers: $num_of_heterozygous_select";
    message "Marker comparison times: $total_diff_cal_times"; 
    hr;
}

main();
