#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;

# ------------------------------------------------------------
# summarize_qtl_intervals_from_logs.pl
#
# Input:
#   qtl_trait_report_<TRAIT>_*.txt (per-trait sink output from qtl-mapping-May2026.R)
#   Legacy: qtl_log_<TRAIT>_*.log in the same directory is still recognized.
#
# Output:
#   one row per QTL, using intervals printed in the log under:
#       "Intervals by chromosome (from scanone output):"
#
#   plus refined peak position, additive/dominance effects,
#   QTL %var and Pvalue(F), and F1parent %var from:
#       "refineqtl + refit"
#
# Usage (from qtl_mapping repo root):
#   perl summarize_qtl_intervals_from_logs.pl \
#     --log_dir results/trait_qtl_reports \
#     --interval lod_1.5 \
#     --out qtl_lod15_summary.tsv
#
# Valid intervals:
#   lod_1.5
#   lod_2
#   bayes
# ------------------------------------------------------------

my %args = @ARGV;

my $log_dir  = $args{"--log_dir"} // "results/trait_qtl_reports";
my $out_file = $args{"--out"}     // "qtl_interval_summary.tsv";

my $interval_to_report = $args{"--interval"} // "lod_1.5";

die "Invalid --interval '$interval_to_report'. Use lod_1.5, lod_2, or bayes.\n"
    unless $interval_to_report =~ /^(lod_1\.5|lod_2|bayes)$/;

opendir(my $dh, $log_dir) or die "Cannot open log_dir '$log_dir': $!\n";
my @log_files = sort grep {
    /^(?:qtl_trait_report_.*\.txt|qtl_log_.*\.log)$/
} readdir($dh);
closedir($dh);

die "No qtl_trait_report_*.txt or qtl_log_*.log files found in '$log_dir'\n" unless @log_files;

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
        qtl_pve_percent
        qtl_pvalue_f
        f1parent_pve_percent
        trait_report_file
        notes
    )
), "\n";

foreach my $log_file (@log_files) {
    my $path = "$log_dir/$log_file";

    my @runs = parse_log_runs($path, $log_file);

    foreach my $run (@runs) {
        my $trait = $run->{trait} // infer_trait_from_log_name($log_file);

        my @qtls = @{ $run->{interval_qtls} // [] };

        if (!@qtls) {
            warn "No '$interval_to_report' intervals found for trait=$trait in $log_file\n";
            next;
        }

        for my $i (0 .. $#qtls) {
            my $qtl_num = $i + 1;
            my $qtl_id  = "Q$qtl_num";
            my $q       = $qtls[$i];

            my $chr_clean = $q->{chr};
            $chr_clean =~ s/^chr//;

            my @notes;

            my $refined_pos = "";
            my $additive    = "";
            my $dominance   = "";
            my $direction   = "";
            my $qtl_pve     = "";
            my $qtl_pval_f  = "";
            my $f1_pve      = $run->{f1parent_pve_percent} // "";

            if (exists $run->{refined_qtls}{$qtl_id}) {
                my $r = $run->{refined_qtls}{$qtl_id};

                $refined_pos = $r->{pos}       // "";
                $additive    = $r->{additive}  // "";
                $dominance   = $r->{dominance} // "";
                $qtl_pve     = $r->{pve}       // "";
                $qtl_pval_f  = $r->{pvalue_f}  // "";

                if ($additive ne "") {
                    $direction = $additive > 0 ? "B_increases_trait" :
                                 $additive < 0 ? "B_decreases_trait" :
                                                 "no_additive_effect";
                }

                if (defined $r->{chr} && $r->{chr} ne $chr_clean) {
                    push @notes, "WARNING_chr_mismatch_interval_chr=$chr_clean refined_chr=$r->{chr}";
                }
            } else {
                push @notes, "No refined log info found for $qtl_id";
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
                $qtl_pve,
                $qtl_pval_f,
                $f1_pve,
                $log_file,    # basename of trait report (.txt or legacy .log)
                join(";", @notes)
            ), "\n";
        }
    }
}

close $OUT;

print "Wrote: $out_file\n";


# ------------------------------------------------------------
# parse_log_runs
#
# Handles either:
#   1. logs with multiple "---- Trait run ----" sections, or
#   2. a single-trait log without explicit per-trait splitting.
# ------------------------------------------------------------
sub parse_log_runs {
    my ($path, $log_file) = @_;

    open(my $IN, "<", $path) or die "Cannot read log file '$path': $!\n";
    local $/ = undef;
    my $txt = <$IN>;
    close $IN;

    my @runs;

    # Split multi-trait logs at "---- Trait run ----".
    if ($txt =~ /---- Trait run ----/) {
        my @pieces = split /(?=---- Trait run ----)/, $txt;

        foreach my $piece (@pieces) {
            next unless $piece =~ /---- Trait run ----/;
            next unless $piece =~ /trait=([A-Za-z0-9_.-]+)/;

            my $trait = $1;
            my $run = parse_one_run_block($piece, $trait);
            push @runs, $run if $run;
        }
    } else {
        my $trait = infer_trait_from_log_name($log_file);

        if ($txt =~ /trait=([A-Za-z0-9_.-]+)/) {
            $trait = $1;
        } elsif ($txt =~ /scanone complete for trait=([A-Za-z0-9_.-]+)/) {
            $trait = $1;
        }

        my $run = parse_one_run_block($txt, $trait);
        push @runs, $run if $run;
    }

    return @runs;
}


# ------------------------------------------------------------
# parse_one_run_block
# ------------------------------------------------------------
sub parse_one_run_block {
    my ($txt, $trait) = @_;

    my %run = (
        trait                  => $trait,
        interval_qtls          => [],
        refined_qtls           => {},
        f1parent_pve_percent   => "",
    );

    my @interval_qtls = parse_intervals_from_log($txt);
    $run{interval_qtls} = \@interval_qtls;

    my %refined = parse_refined_qtls_from_log($txt);
    my %effects = parse_effect_estimates_from_log($txt);

    foreach my $qtl_id (keys %refined) {
        my $chr      = $refined{$qtl_id}{chr};
        my $pos      = $refined{$qtl_id}{pos};
        my $qtl_name = $refined{$qtl_id}{name};

        $refined{$qtl_id}{additive}  = "";
        $refined{$qtl_id}{dominance} = "";

        if (defined $qtl_name && $qtl_name ne "") {
            $refined{$qtl_id}{additive}  = $effects{$qtl_name . "a"} // "";
            $refined{$qtl_id}{dominance} = $effects{$qtl_name . "d"} // "";
        }

        # Fallback: effect labels usually use one decimal place.
        if ($refined{$qtl_id}{additive} eq "" && defined $chr && defined $pos) {
            my $short_pos = sprintf("%.1f", $pos);
            $refined{$qtl_id}{additive}  = $effects{"$chr\@$short_pos" . "a"} // "";
            $refined{$qtl_id}{dominance} = $effects{"$chr\@$short_pos" . "d"} // "";
        }
    }

    my %drop = parse_refined_drop_one_from_log($txt);

    foreach my $qtl_id (keys %refined) {
        my $qtl_name = $refined{$qtl_id}{name};

        if (defined $qtl_name && exists $drop{$qtl_name}) {
            $refined{$qtl_id}{pve}      = $drop{$qtl_name}{pve};
            $refined{$qtl_id}{pvalue_f} = $drop{$qtl_name}{pvalue_f};
        }
    }

    if (exists $drop{F1parent}) {
        $run{f1parent_pve_percent} = $drop{F1parent}{pve};
    }

    $run{refined_qtls} = \%refined;

    return \%run;
}


# ------------------------------------------------------------
# parse_intervals_from_log
#
# Pulls rows under:
#   Intervals by chromosome (from scanone output):
#
# Keeps only the interval requested globally by $interval_to_report.
# Groups by chr. Within each chr:
#   start = lowest position
#   end   = highest position
#   peak  = highest LOD
# ------------------------------------------------------------
sub parse_intervals_from_log {
    my ($txt) = @_;

    my @blocks = ($txt =~ /Intervals by chromosome \(from scanone output\):\s*\n(.*?)(?=\n\[.*?\]\s+Wrote intervals table|\n----|\z)/sg);

    return () unless @blocks;

    # Use the last interval block in this run.
    my $block = $blocks[-1];

    my %by_chr;
    my @chr_order;
    my %seen_chr;

    foreach my $line (split /\n/, $block) {
        $line =~ s/^\s+|\s+$//g;
        next if $line eq "";
        next if $line =~ /^chr\s+interval\s+marker\s+pos\s+lod$/;
        next unless $line =~ /^\d+\s+/;

        my @f = split /\s+/, $line;

        # Expected:
        # row_index chr interval marker pos lod
        # 1 chr2 lod_1.5 Chrom2-18510654 38.49010 2.445721
        next unless @f >= 6;

        my $chr      = $f[1];
        my $interval = $f[2];
        my $marker   = $f[3];
        my $pos      = $f[4];
        my $lod      = $f[5];

        next unless $interval eq $interval_to_report;

        if (!$seen_chr{$chr}++) {
            push @chr_order, $chr;
        }

        push @{ $by_chr{$chr} }, {
            marker => $marker,
            pos    => $pos,
            lod    => $lod,
        };
    }

    my @qtls;

    foreach my $chr (@chr_order) {
        my @rows = @{ $by_chr{$chr} };

        if (@rows < 3) {
            warn "Expected at least 3 '$interval_to_report' rows for trait interval $chr; found " . scalar(@rows) . "\n";
            next;
        }

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
# parse_refined_qtls_from_log
#
# Pulls QTL names/chrs/positions from:
#   Refined QTL object:
#
# Expected row:
#   Q1 2@46.9 2 46.853 3
# ------------------------------------------------------------
sub parse_refined_qtls_from_log {
    my ($txt) = @_;

    my %refined;

    my @blocks = ($txt =~ /Refined QTL object:\s*\n(.*?)(?=\n\[.*?\]\s+Refined model PVE|\n----|\z)/sg);
    return %refined unless @blocks;

    my $block = $blocks[-1];

    foreach my $line (split /\n/, $block) {
        $line =~ s/^\s+|\s+$//g;
        next unless $line =~ /^Q\d+\s+/;

        my @f = split /\s+/, $line;
        next unless @f >= 4;

        my $qtl_id = $f[0];
        my $name   = $f[1];
        my $chr    = $f[2];
        my $pos    = $f[3];

        $refined{$qtl_id} = {
            name => $name,
            chr  => $chr,
            pos  => $pos,
        };
    }

    return %refined;
}


# ------------------------------------------------------------
# parse_effect_estimates_from_log
#
# Pulls additive/dominance effects from:
#   Effect estimates:
#
# Expected:
#   Intercept F1parent 2@46.9a 2@46.9d ...
#   14.04     0.62     11.91   0.36 ...
# ------------------------------------------------------------
sub parse_effect_estimates_from_log {
    my ($txt) = @_;

    my %effects;

    my @blocks = ($txt =~ /Effect estimates:\s*\n\s*(.+?)\n\s*(.+?)(?=\n\n|\n----|\z)/sg);
    return %effects unless @blocks >= 2;

    my $values_line = pop @blocks;
    my $names_line  = pop @blocks;

    $names_line  =~ s/^\s+|\s+$//g;
    $values_line =~ s/^\s+|\s+$//g;

    my @names  = split /\s+/, $names_line;
    my @values = split /\s+/, $values_line;

    for my $i (0 .. $#names) {
        $effects{$names[$i]} = $values[$i] if defined $values[$i];
    }

    return %effects;
}


# ------------------------------------------------------------
# parse_refined_drop_one_from_log
#
# Pulls the Drop-one table that occurs after:
#   "refineqtl + refit"
#
# Expected rows:
#   2@46.9    2    10982.01 6.799221 12.268666 ... Pvalue(F)
#   F1parent  1     8042.31 5.081955  8.984548 ... Pvalue(F)
#
# Returns:
#   $drop{"2@46.9"}{pve}
#   $drop{"2@46.9"}{pvalue_f}
#   $drop{"F1parent"}{pve}
# ------------------------------------------------------------
sub parse_refined_drop_one_from_log {
    my ($txt) = @_;

    my %drop;

    # First restrict to the refineqtl + refit section.
    my @sections = ($txt =~ /---- refineqtl \+ refit ----(.*?)(?=\n----|\z)/sg);
    return %drop unless @sections;

    my $section = $sections[-1];

    my @blocks = ($section =~ /Drop-one table:\s*\n(.*?)(?=\nattr\(|\n\[.*?\]\s+Effect estimates:|\z)/sg);
    return %drop unless @blocks;

    my $block = $blocks[-1];

    foreach my $line (split /\n/, $block) {
        $line =~ s/^\s+|\s+$//g;
        next if $line eq "";
        next if $line =~ /^df\s+Type\s+III\s+SS/;
        next if $line =~ /^attr\(/;

        my @f = split /\s+/, $line;

        # Expected after split:
        # marker df Type III SS LOD %var F value Pvalue(Chi2) Pvalue(F)
        #
        # Because "Type III SS", "F value", and "Pvalue(...)" print as
        # separated text in the header, but the data rows are simple:
        # 2@46.9 2 10982.01 6.799221 12.268666 16.44761 1.58e-07 2.55e-07
        next unless @f >= 8;

        my $term     = $f[0];
        my $pve      = $f[4];
        my $pvalue_f = $f[7];

        next unless $term =~ /^(?:F1parent|[\d.]+\@[-+]?\d+(?:\.\d+)?)$/;

        $drop{$term} = {
            pve      => $pve,
            pvalue_f => $pvalue_f,
        };
    }

    return %drop;
}


# ------------------------------------------------------------
# infer_trait_from_log_name
# ------------------------------------------------------------
sub infer_trait_from_log_name {
    my ($log_file) = @_;

    my $trait = $log_file;
    if ($trait =~ /^qtl_trait_report_(.+)_\d{8}_\d{6}\.txt$/) {
        return $1;
    }
    $trait =~ s/^qtl_log_//;
    $trait =~ s/_\d{8}_\d{6}\.log$//;
    $trait =~ s/\.log$//;

    return $trait;
}
