#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use List::Util qw(max);

# ------------------------------------------------------------
# summarize_qtl_intervals.pl
#
# Input:
#   *_peak_intervals.tsv
#   qtl_log_<TRAIT>_*.log
#
# Output:
#   one row per QTL, using lod_1.5 intervals, plus refined
#   peak position and effect estimates from the matching log file.
#
# Usage:
#   perl summarize_qtl_intervals.pl \
#     --interval_dir qtl_mapping \
#     --log_dir qtl_mapping \
#     --interval lod_1.5|lod_2|bayes
#     --out qtl_lod15_summary.tsv
# ------------------------------------------------------------

my %args = @ARGV;

my $interval_dir = $args{"--interval_dir"} // ".";
my $log_dir      = $args{"--log_dir"}      // ".";
my $out_file     = $args{"--out"}          // "qtl_lod15_summary.tsv";

my $interval_to_report = $args{"--interval"} // "lod_1.5";

die "Invalid --interval '$interval_to_report'. Use lod_1.5, lod_2, or bayes.\n"
    unless $interval_to_report =~ /^(lod_1\.5|lod_2|bayes)$/;

opendir(my $dh, $interval_dir) or die "Cannot open interval_dir '$interval_dir': $!\n";
my @interval_files = sort grep { /_peak_intervals.*\.tsv$/ } readdir($dh);
closedir($dh);

open(my $OUT, ">", $out_file) or die "Cannot write '$out_file': $!\n";

print $OUT join("\t",
    qw(
        trait
        qtl_id
        chr
        interval_type
        start_marker
        start_pos
        start_lod
        peak_marker
        peak_pos
        peak_lod
        end_marker
        end_pos
        end_lod
        refined_peak_pos
        additive
        dominance
        B_allele_direction
        log_file
        interval_file
        notes
    )
), "\n";

foreach my $file (@interval_files) {
    my $path = "$interval_dir/$file";

    my $trait = $file;
    $trait =~ s/_peak_intervals.*\.tsv$//;

    my $log_file = find_matching_log($log_dir, $trait);

    my $log_info = {};
    if ($log_file) {
        $log_info = parse_log_file("$log_dir/$log_file");
    } else {
        warn "No log file found for trait=$trait\n";
    }

    my @qtls = parse_intervals($path, $interval_to_report);

    for my $i (0 .. $#qtls) {
        my $qtl_num = $i + 1;
        my $qtl_id  = "Q$qtl_num";
        my $q       = $qtls[$i];

        my $chr_clean = $q->{chr};
        $chr_clean =~ s/^chr//;

        my $refined_pos = "";
        my $additive    = "";
        my $dominance   = "";
        my $direction   = "";
        my @notes;

        if (exists $log_info->{$qtl_id}) {
            $refined_pos = $log_info->{$qtl_id}->{pos} // "";
            $additive    = $log_info->{$qtl_id}->{additive} // "";
            $dominance   = $log_info->{$qtl_id}->{dominance} // "";

            if ($additive ne "") {
                $direction = $additive > 0 ? "B_increases_trait" :
                             $additive < 0 ? "B_decreases_trait" :
                                             "no_additive_effect";
            }
        } else {
            push @notes, "No refined log info found for $qtl_id";
        }

        # Sanity check: QTL order should generally match between interval file
        # and refined QTL object. This warns if chromosome order looks mismatched.
        if (exists $log_info->{$qtl_id} && defined $log_info->{$qtl_id}->{chr}) {
            my $log_chr = $log_info->{$qtl_id}->{chr};
            if ($log_chr ne $chr_clean) {
                push @notes, "WARNING_chr_mismatch_interval_chr=$chr_clean log_chr=$log_chr";
            }
        }

        print $OUT join("\t",
            $trait,
            $qtl_id,
            $q->{chr},
            $interval_to_report,

            $q->{start_marker},
            $q->{start_pos},
            $q->{start_lod},

            $q->{peak_marker},
            $q->{peak_pos},
            $q->{peak_lod},

            $q->{end_marker},
            $q->{end_pos},
            $q->{end_lod},

            $refined_pos,
            $additive,
            $dominance,
            $direction,
            $log_file // "",
            $file,
            join(";", @notes)
        ), "\n";
    }
}

close $OUT;

print "Wrote: $out_file\n";


# ------------------------------------------------------------
# Parse one *_peak_intervals.tsv file.
# Keeps only interval == lod_1.5.
# Assumes each QTL/chromosome has three lod_1.5 rows:
#   start, peak, end
# Peak is identified as row with maximum LOD.
# ------------------------------------------------------------
sub parse_intervals {
    my ($path, $interval_to_report) = @_;

    open(my $IN, "<", $path) or die "Cannot read '$path': $!\n";

    my $header = <$IN>;
    chomp $header;
    my @cols = split /\t/, $header;
    my %idx;
    for my $i (0 .. $#cols) {
        $idx{$cols[$i]} = $i;
    }

    for my $needed (qw(chr interval marker pos lod)) {
        die "Missing required column '$needed' in $path\n"
            unless exists $idx{$needed};
    }

    my %by_chr;
    my @chr_order;
    my %seen;

    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        my @f = split /\t/, $line;

        next unless $f[$idx{interval}] eq $interval_to_report;

        my $chr = $f[$idx{chr}];

        if (!$seen{$chr}++) {
            push @chr_order, $chr;
        }

        push @{ $by_chr{$chr} }, {
            marker => $f[$idx{marker}],
            pos    => $f[$idx{pos}],
            lod    => $f[$idx{lod}],
        };
    }

    close $IN;

    my @qtls;

    foreach my $chr (@chr_order) {
        my @rows = @{ $by_chr{$chr} };

        if (@rows < 3) {
            warn "Expected at least 3 '$interval_to_report' rows for $path $chr; found " . scalar(@rows) . "\n";
            next;
        }

        # Sort by map position so start/end are position-based.
        @rows = sort { $a->{pos} <=> $b->{pos} } @rows;

        my $start = $rows[0];
        my $end   = $rows[-1];

        my $peak = $rows[0];
        foreach my $r (@rows) {
            $peak = $r if $r->{lod} > $peak->{lod};
        }

        push @qtls, {
            chr          => $chr,

            start_marker => $start->{marker},
            start_pos    => $start->{pos},
            start_lod    => $start->{lod},

            peak_marker  => $peak->{marker},
            peak_pos     => $peak->{pos},
            peak_lod     => $peak->{lod},

            end_marker   => $end->{marker},
            end_pos      => $end->{pos},
            end_lod      => $end->{lod},
        };
    }

    return @qtls;
}

# ------------------------------------------------------------
# Find the most recent matching qtl_log_<TRAIT>_*.log file.
# If you have multiple logs per trait, this chooses the newest
# by file modification time.
# ------------------------------------------------------------
sub find_matching_log {
    my ($log_dir, $trait) = @_;

    opendir(my $dh, $log_dir) or die "Cannot open log_dir '$log_dir': $!\n";
    my @logs = grep { /^qtl_log_\Q$trait\E_.*\.log$/ } readdir($dh);
    closedir($dh);

    return undef unless @logs;

    @logs = sort {
        (stat("$log_dir/$b"))[9] <=> (stat("$log_dir/$a"))[9]
    } @logs;

    return $logs[0];
}


# ------------------------------------------------------------
# Parse one QTL log.
#
# Pulls from the final/refined section:
#   Refined QTL object:
#      Q1 2@46.9   2 46.853 3
#      Q2 1@94.5   1 94.500 3
#
# Then pulls additive/dominance terms from:
#   Effect estimates:
#      Intercept F1parent 2@46.9a 2@46.9d ...
#      values...
#
# Returns:
#   $info->{Q1}->{chr}
#   $info->{Q1}->{pos}
#   $info->{Q1}->{additive}
#   $info->{Q1}->{dominance}
# ------------------------------------------------------------
sub parse_log_file {
    my ($path) = @_;

    open(my $IN, "<", $path) or die "Cannot read log file '$path': $!\n";
    local $/ = undef;
    my $txt = <$IN>;
    close $IN;

    my %info;

    # Use the last Refined QTL object in the file, in case logs contain multiple runs.
    my @refined_blocks = ($txt =~ /\[.*?\]\s+Refined QTL object:\s*\n(.*?)(?=\n\[.*?\] Refined model PVE|\n----|\z)/sg);

    if (@refined_blocks) {
        my $block = $refined_blocks[-1];

        foreach my $line (split /\n/, $block) {
            $line =~ s/^\s+|\s+$//g;
            next unless $line =~ /^Q\d+\s+/;

            my @f = split /\s+/, $line;

            # Expected:
            # Q1 2@46.9 2 46.853 3
            my $qtl_id = $f[0];
            my $chr    = $f[2];
            my $pos    = $f[3];

            $info{$qtl_id}->{chr} = $chr;
            $info{$qtl_id}->{pos} = $pos;
        }
    }

    # Use the last Effect estimates section in the file.
    my @effect_blocks = ($txt =~ /Effect estimates:\s*\n\s*(.+?)\n\s*(.+?)(?=\n\n|\n----|\z)/sg);

    if (@effect_blocks >= 2) {
        my $values_line = pop @effect_blocks;
        my $names_line  = pop @effect_blocks;

        $names_line  =~ s/^\s+|\s+$//g;
        $values_line =~ s/^\s+|\s+$//g;

        my @names  = split /\s+/, $names_line;
        my @values = split /\s+/, $values_line;

        my %effects;
        for my $i (0 .. $#names) {
            $effects{$names[$i]} = $values[$i] if defined $values[$i];
        }

        foreach my $qtl_id (keys %info) {
            my $chr = $info{$qtl_id}->{chr};
            my $pos = $info{$qtl_id}->{pos};

            next unless defined $chr && defined $pos;

            my $short_pos = sprintf("%.1f", $pos);

            # Effect labels usually look like:
            #   2@46.9a
            #   2@46.9d
            #   1@94.5a
            #   1@94.5d
            #
            # Because R/qtl rounds names in the printed table, use one decimal.
            my $a_key = "$chr\@$short_pos" . "a";
            my $d_key = "$chr\@$short_pos" . "d";

            $info{$qtl_id}->{additive}  = $effects{$a_key} // "";
            $info{$qtl_id}->{dominance} = $effects{$d_key} // "";
        }
    }

    return \%info;
}