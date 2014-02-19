#!/usr/bin/env python
"""
precision_recall_curves.py

Writes PDF plot of precision vs. recall of a spliced aligner as each of a set
of assorted parameters are varied given a BED file T with true introns
(typically from a simulator like Flux) and a BED file Y with introns
retrieved by the aligner. Y stores reported coverage in the score field and
other parameters in a semicolon-separated list in the name field. For example,
the name field of a line from Y should look like

JUNC0001;anchor_significance=8;maximum_match_rate=.89;
unique_displacement_count=4

. Also plots and outputs values of parameters at which (precision + recall)
is maximized along the curves.

A line in an input BED file characterizes an intron by decking it with blocks.
A single intron defines two junctions, one on either side of it. A true
positive is an intron recovered by the aligner, a false positive is an intron
that appears in Y but not in T, and a false negative is an intron that appears
in T but not in Y. Precision is
(# true positives) / (# true positives + # false positives) . Recall is
(# true positives) / (# true positives + # false negatives) .

Check out http://www.huyng.com/posts/sane-color-scheme-for-matplotlib/ to
make matplotlib plots prettier.

Output (written to stdout)
----------------------------
Values of parameters at which precision + recall is maximized along curves
in tab-delimited format:
1. Parameter name (e.g., coverage)
2. Parameter p's value (argmax (precision(p) + recall(p)))
3. Precision at max
4. Recall at max
"""

import sys
import matplotlib.pyplot as plt
import string

def introns_from_bed_stream(bed_stream, read_stats=False):
    """ Converts BED to dictionary that maps RNAMES to sorted lists of introns.

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions.
        read_stats: True iff statistics should be read. Statistics
        should be in a semicolon-separated list in the name field --
        (e.g., JUNC0001;anchor_significance=8;maximum_match_rate=.89;
        unique_displacement_count=4). Coverage is taken to be in the score
        field. read_stats=True assumes there is only one intron (two blocks)
        per line. In general read_statistics=False when reading a Flux BED.

        Return value: if read_stats=False: a dictionary. Each key is an RNAME,
            typically a chromosome, and its corresponding value is a sorted
            list of unique introns. An intron is specified by a tuple
            (pos, end_pos) The primary sort key is pos, and the
            secondary sort key is end_pos.
            if readstats=True: a tuple with two dictionaries: 1) each
            key is an RNAME, typically a chromosome, and its corresponding
            value is a sorted list of unique introns. An intron is specified by
            a tuple (pos, end_pos, stats), where pos is the start position
            (inclusive) and end_pos is the end position (exclusive). stats is a
            dictionary whose keys are statistic names and whose values are
            their corresponding values. 2) each key is a statistic, and its
            corresponding value is a tuple (discrete, min, max), where min/max
            are the min/max values of the statistic of all introns processed,
            and discrete is True if the statistic is discrete; otherwise, it's
            False.
    """
    # For separating stat names into words
    underscore_to_space = string.maketrans('_', ' ')
    introns = {}
    if read_stats:
        # Store maxes and mins of stats as well as whether stats are discrete
        stat_ranges = {}
        stat_discrete = {}
    for line in bed_stream:
        tokens = line.rstrip().split('\t')
        if len(tokens) != 12:
            continue
        chrom = tokens[0]
        chrom_start = int(tokens[1])
        chrom_end = int(tokens[2])
        name = tokens[3]
        if read_stats:
            stats = {}
            for stat in name.rstrip().split(';')[1:]:
                stat_tokens = stat.split('=')
                stat_name = stat_tokens[0].translate(underscore_to_space)
                stats[stat_name] = float(stat_tokens[1])
                int_stat = int(stats[stat_name])
                if stat_name not in stat_discrete:
                    stat_discrete[stat_name] = True
                if int_stat == stats[stat_name]:
                    # Convert to integer
                    stats[stat_name] = int_stat
                else:
                    stat_discrete[stat_name] = False
                if stat_name not in stat_ranges:
                    stat_ranges[stat_name] = [stats[stat_name],
                                                stats[stat_name]]
                else:
                    stat_ranges[stat_name][0] = min(stat_ranges[stat_name][0],
                                                        stats[stat_name])
                    stat_ranges[stat_name][1] = max(stat_ranges[stat_name][1],
                                                        stats[stat_name])
            stats['coverage'] = int(tokens[4])
            if 'coverage' not in stat_ranges:
                stat_ranges['coverage'] = [stats['coverage'],
                                            stats['coverage']]
            else:
                stat_ranges['coverage'][0] = min(stat_ranges['coverage'][0],
                                                    stats['coverage'])
                stat_ranges['coverage'][1] = max(stat_ranges['coverage'][1],
                                                    stats['coverage'])
        if chrom not in introns:
            # Use set so there are no duplicate introns when reading Flux data
            introns[chrom] = ([] if read_stats else set())
        block_sizes = tokens[10].split(',')
        block_starts = tokens[11].split(',')
        # Handle trailing commas
        try:
            int(block_sizes[-1])
        except ValueError:
            block_sizes = block_sizes[:-1]
        try:
            int(block_starts[-1])
        except ValueError:
            block_starts = block_starts[:-1]
        block_count = len(block_sizes)
        if block_count < 2:
            # No introns
            continue
        assert block_count == len(block_starts)
        junctions = []
        # First block characterizes junction on left side of intron
        junctions.append(chrom_start + int(block_starts[0]) 
                                + int(block_sizes[0]))
        for i in xrange(1, block_count - 1):
            # Any intervening blocks characterize two junctions
            intron_start = chrom_start + int(block_starts[i])
            junctions.append(intron_start)
            junctions.append(intron_start + int(block_sizes[i]))
        # Final block characterizes junction on right side of intron
        junctions.append(chrom_start + int(block_starts[-1]))
        # Now associate junctions to call introns
        for i in xrange(len(junctions) - 1):
            if read_stats:
                introns[chrom].append((junctions[i], junctions[i+1], stats))
            else:
                introns[chrom].add((junctions[i], junctions[i+1]))
    if read_stats:
        for chrom in introns:
            introns[chrom].sort()
        to_return = {}
        for stat in stat_ranges:
            if stat == 'coverage': continue
            to_return[stat] = tuple([stat_discrete[stat]]) \
                                + tuple(stat_ranges[stat])
        to_return['coverage'] = tuple([True]) \
                                  + tuple(stat_ranges['coverage'])
        return introns, to_return
    else:
        for chrom in introns:
            introns[chrom] = sorted(list(introns[chrom]))
        return introns

def filtered_introns(introns, stat_coefficients, threshold):
    """ Removes introns for which superposition of stats is below a threshold.

        For a given intron, let S = sum_i a_i s_i, where a_i is the coefficient
        of statistic s_i as specified in stat_coefficients. If S falls below
        threshold, the intron is filtered out.

        introns: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a list of
            introns. An intron is specified by a tuple (pos, end_pos, stats),
            where pos is the start position (inclusive) and end_pos is the end
            position (exclusive). stats is a dictionary whose keys are
            statistic names and whose values are their corresponding values.
        stat_coefficients: a dictionary whose keys are statistic names (here,
            underscores should be replaced with spaces) and whose values are
            coefficients. The user may find it convenient to arrange that the
            sum of the a_i is 1.
        threshold: the value of S below which an intron is filtered out.

        Return value: a dictionary just like the input "introns," except every
            intron has statistic S >= threshold, and no statistics are
            tacked onto output tuples. So each intron is just a tuple
            (pos, end_pos).
    """
    filtered_introns = {}
    for chrom, chrom_introns in introns.items():
        if chrom not in filtered_introns:
            filtered_introns[chrom] = []
        for intron in chrom_introns:
            if sum([stat_coefficients[stat] * intron[2][stat]
                    for stat in stat_coefficients]) >= threshold:
                filtered_introns[chrom].append(intron[:-1])
    return filtered_introns

def information_retrieval_metrics(true_introns, retrieved_introns):
    """ Computes IR stats given sorted lists of true and retrieved introns.

        true_introns: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a sorted list of
            unique true introns, each specified by a tuple (pos, end_pos)
            indicating the intron spans [pos, end_pos). The primary sort key is
            pos, and the secondary sort key is end_pos.
        retrieved_introns: a dictionary. Each key is an RNAME, typically a
            chromosome, and its corresponding value is a sorted list of
            unique retrieved introns, each specified by a tuple (pos, end_pos)
            indicating the intron spans [pos, end_pos). The primary sort key is
            pos, and the secondary sort key is end_pos.

        Return value: tuple (recall, precision, true positive count,
                              false positive count, false negative count).
    """
    true_positive_count = 0
    total_retrieved_introns = 0
    for chrom in retrieved_introns:
        if chrom in true_introns:
            true_introns_count = len(true_introns[chrom])
        else:
            true_introns_count = 0
        retrieved_introns_count = len(retrieved_introns[chrom])
        total_retrieved_introns += retrieved_introns_count
        i, j = 0, 0
        while True:
            if j == true_introns_count or i == retrieved_introns_count:
                break
            if retrieved_introns[chrom][i] == true_introns[chrom][j]:
                true_positive_count += 1
                i += 1
                j += 1
            elif retrieved_introns[chrom][i][0] == true_introns[chrom][j][0]:
                if retrieved_introns[chrom][i][1] > true_introns[chrom][j][1]:
                    j += 1
                else:
                    i += 1
            elif retrieved_introns[chrom][i][0] > true_introns[chrom][j][0]:
                j += 1
            else:
                i += 1
    total_true_introns = 0
    for chrom in true_introns:
        total_true_introns += len(true_introns[chrom])
    false_positive_count = total_retrieved_introns - true_positive_count
    false_negative_count = total_true_introns - true_positive_count
    if false_negative_count == 0:
        recall = 1
    else:
        recall = float(true_positive_count) / total_true_introns
    if false_positive_count == 0:
        precision = 1
    else:
        precision = float(true_positive_count) / total_retrieved_introns
    return (recall, precision, true_positive_count, false_positive_count,
                false_negative_count)

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-t', '--true-introns-bed', type=str, required=True, 
        help='Full path of BED file containing true introns. This should be '
             'a Flux BED.')
    parser.add_argument('-r', '--retrieved-introns-bed', type=str,
        required=True, 
        help='Full path of BED file containing introns retrieved by aligner')
    parser.add_argument('-d', '--plot-density', type=int, required=False,
        default=50,
        help='Number of points to compute per curve; the larger this number, '
             'the better the interpolations. If the curve is discrete, all '
             'possible values are plotted, so this number is ignored.')
    parser.add_argument('-o', '--output-filename', type=str, required=False,
        default='precision_recall_curve.pdf',
        help='Filename for output plot. Extension specifies format. pdf and '
             'png are two good choices. See matplotlib documentation for '
             'more.')
    args = parser.parse_args(sys.argv[1:])
    with open(args.true_introns_bed) as true_introns_bed_stream:
        print >>sys.stderr, 'Reading true introns BED....'
        true_introns = introns_from_bed_stream(true_introns_bed_stream)
    with open(args.retrieved_introns_bed) as retrieved_introns_bed_stream:
        print >>sys.stderr, 'Reading retrieved introns BED....'
        retrieved_introns, retrieved_intron_stats \
            = introns_from_bed_stream(retrieved_introns_bed_stream,
                                        read_stats=True)
    print >>sys.stderr, 'Computing precisions and recalls....'
    points_and_thresholds = {}
    import itertools
    plot_series = []
    figure = plt.figure()
    colors = itertools.cycle(['r', 'g', 'b', 'k', 'c', 'm', 'y'])
    min_recall = 1
    max_recall = 0
    min_precision = 1
    max_precision = 0
    for stat, discrete_and_range in retrieved_intron_stats.items():
        discrete = discrete_and_range[0]
        stat_min, stat_max = discrete_and_range[1:]
        if discrete:
            thresholds = range(stat_min, stat_max+1)
        else:
            interval = float(stat_max - stat_min) / (args.plot_density - 1)
            thresholds = [stat_min + interval*i for i in 
                                    xrange(args.plot_density)]
        stat_coefficients = {}
        for stat_name in retrieved_intron_stats.keys():
            stat_coefficients[stat_name] = 0
        stat_coefficients[stat] = 1
        points_and_thresholds[stat] = [information_retrieval_metrics(
                                                true_introns, 
                                                filtered_introns(
                                                        retrieved_introns, 
                                                        stat_coefficients, 
                                                        threshold
                                                    )       
                                            ) + (threshold,)
                                        for threshold in thresholds]
        recalls = [point_and_threshold[0] for point_and_threshold in
                        points_and_thresholds[stat]]
        precisions = [point_and_threshold[1] for point_and_threshold in
                        points_and_thresholds[stat]]
        min_recall = min(min_recall, min(recalls))
        max_recall = max(max_recall, max(recalls))
        min_precision = min(min_precision, min(precisions))
        max_precision = max(max_precision, max(precisions))
        current_color=next(colors)
        plot_series.append(
                plt.scatter(recalls, precisions, c=current_color, s=75)
            )
        plt.plot(recalls, precisions, c=current_color)
    plt.legend(plot_series, 
                retrieved_intron_stats.keys(),
                scatterpoints=1,
                loc='lower right',
                fontsize=16)
    plt.xlabel('Recall', fontsize=20)
    plt.ylabel('Precision', fontsize=20)
    padding = .015
    x_padding = (max_recall - min_recall) * padding
    y_padding = (max_precision - min_precision) * padding
    plt.xlim(min_recall - x_padding, max_recall + x_padding)
    plt.ylim(min_precision - y_padding, max_precision + y_padding)
    figure.suptitle('Precision and recall for intron calls satisfying '
                    'various parameter thresholds',
        fontsize=20)
    figure.savefig(args.output_filename)
    # Print best values
    for stat in points_and_thresholds:
        highest_precision_recall_sum = None
        best_point_and_threshold = None
        for point_and_threshold in points_and_thresholds[stat]:
            current_sum = sum(point_and_threshold[:2])
            if current_sum > highest_precision_recall_sum:
                highest_precision_recall_sum = current_sum
                best_recall = point_and_threshold[0]
                best_precision = point_and_threshold[1]
                best_point_and_threshold = point_and_threshold
        print '%s\t%s\t%.12f\t%.12f' % (stat, 
                                        str(best_point_and_threshold[-1]),
                                        best_precision, best_recall)
