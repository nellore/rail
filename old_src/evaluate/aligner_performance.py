#!/usr/bin/env python
"""
aligner_performance.py

Reads multiple BED files and draws graphs of sensitivity and specificity as a function of junction coverage.

Input (command line arguments)
________________________________
At least 2 BED files. The first file is the junction reference file. It should be in the same format
as the junctions.bed output from a flux simulation. All other BED files should be in the format output 
by the rail pipeline (specifically including a coverage value for each junction.)


Output
----------------------------
2 graphs for specificity and sensitivity as a function of junction coverage in the reference file.
Each graph contains a curve for each input BED file other than the reference file.
"""

import sys
from bisect import bisect_left
from operator import itemgetter
import matplotlib.pyplot as plt

def junctions_by_coverage(junctions):
    """
        Convert a list of junctions to a dictionary containing a list of junctions for each coverage level.

        junctions: List of junctions, each represented as a tuple of a chromosome name and index.

        Return value: a dictionary. Each key is a coverage level, and corresponding value is a list of 
            junctions with that level of coverage.
    """
    
    # Sort the junctions
    junctions = sorted(junctions)

    print junctions

    # Loop through junctions and split into lists based on coverage level
    coverage_levels = dict()

    curr_coverage = 0
    curr_junction = junctions[0]
    for j in junctions:
        if j == curr_junction:
            curr_coverage += 1
        else:
            if curr_coverage in coverage_levels:
                coverage_levels[curr_coverage].append(curr_junction)
            else:
                coverage_levels[curr_coverage] = [curr_junction]

            curr_junction = j
            curr_coverage = 1

    # Print coverage distribution
    for k,v in coverage_levels.items():
        print '%d junctions with coverage %d' % (len(v), k)

    return coverage_levels


def true_junctions_from_bed(bed_stream):
    """ Parses the BED file and returns an unsorted list of tuples. Each tuple contains the 
            chromosome and index of a junction from the BED file. May return duplicate junctions
            if they are present in the BED file.

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions. The BED file should be in the same format as the flux
            output junctions file, with a line corresponding to each read.

        Return value: a list of tuples. Each tuple contains a chromosome name and
            index, pointing to a junction location from the BED file, and the coverage
            for that junction.
    """

    print 'Reading true junctions'

    # Reduce the BED file to a list of tuples.
    # Each tuple contains a chromosome and index within that chromosome, corresponding to a junction in a read.

    reads = 0
    spliced = 0

    read_junctions = []
    for line in bed_stream:
        row = line.rstrip().split('\t')

        reads += 1

        if len(row) != 12 or int(row[-3]) == 1:
            continue

        spliced += 1

        segment_lens = [int(x.rstrip()) for x in row[-2].split(',')]
        segment_offsets = [int(x.rstrip()) for x in row[-1].split(',')]

        for i in xrange(len(segment_lens)-1):
            start = int(row[1]) + segment_lens[i] + segment_offsets[i]
            end = int(row[1]) + segment_offsets[i+1]
            
            read_junctions.append((row[0], start, end))

            #read_junctions.append((row[0], index))

            #index = int(row[1]) + segment_offsets[i+1]
            #read_junctions.append((row[0], index))

    read_junctions = sorted(read_junctions)

    #print '%d reads' % reads
    #print '%d spliced' % spliced
    #print '%d junctions' % len(read_junctions)

    junctions = []
    curr_coverage = 0
    curr_junction = read_junctions[0]
    for j in read_junctions:
        if j == curr_junction:
            curr_coverage += 1
        else:
            junctions.append((curr_junction[0], curr_junction[1], curr_junction[2], curr_coverage))

            curr_junction = j
            curr_coverage = 1


    #print '%d unique junctions' % len(junctions)

    return junctions

def retrieved_junctions_from_bed(bed_stream):
    """ Parses the BED file and returns an unsorted list of tuples. Each tuple contains the 
            chromosome and index of a junction from the BED file. May return duplicate junctions
            if they are present in the BED file.

        bed_stream: input stream containing lines of a BED file characterizing
            splice junctions. The BED file should be in the same format as the
            rail output, with a line for each junction.

        Return value: a list of tuples. Each tuple contains a chromosome name and
            index, pointing to a junction location from the BED file, and the coverage
            for that junction.
    """

    print 'Reading retrieved junctions'

    # Reduce the BED file to a list of tuples.
    # Each tuple contains a chromosome and index within that chromosome, corresponding to a junction in a read.
    read_junctions = []
    for line in bed_stream:
        row = line.rstrip().split('\t')
        if len(row) != 12:
            continue

        segment_lens = [int(x.rstrip()) for x in row[-2].split(',')]
        segment_offsets = [int(x.rstrip()) for x in row[-1].split(',')]


        for i in xrange(len(segment_lens)-1):
            start = int(row[1]) + segment_lens[i] + segment_offsets[i]
            end = int(row[1]) + segment_offsets[i+1]

            read_junctions.append((row[0], start, end, int(row[4])))
            #read_junctions.append((row[0], index, int(row[4])))

            #index = int(row[1]) + segment_offsets[i+1]
            #read_junctions.append((row[0], index))

    read_junctions = sorted(read_junctions)

    return read_junctions

    #junctions = []
    #curr_coverage = 0
    #curr_junction = read_junctions[0]
    #for j in read_junctions:
    #    if j == curr_junction:
    #        curr_coverage += 1
    #    else:
    #        junctions.append((curr_junction[0], curr_junction[1], curr_coverage))
    #
    #        curr_junction = j
    #        curr_coverage = 1
    #
    #return junctions


def compare_junctions(true_junctions, retrieved_junctions):
    """ Compare the true and retrieved junctions and return a list, separated by true junction
            coverage, of true positives, false positives, and false negatives.

        true_junctions: Dictionary containing a list of junction tuples as values, with coverage
            as the key.

        retrieved_junctions: List of tuples corresponding to 

    """

    print 'Comparing junctions'

    # Find maximum coverage values
    true_max = max(true_junctions, key=itemgetter(3))[3]
    retrieved_max = max(retrieved_junctions, key=itemgetter(3))[3]

    print 'Max true coverage = %d' % true_max
    print 'Max retrieved coverage = %d' % retrieved_max

    # For recall, count the true positives and total true for true coverage levels
    recall_tp = [0]*(true_max+1)
    recall_t = [0]*(true_max+1)

    # For precision, count the true posiives and total positives for retrieved coverage levels
    precision_tp = [0]*(retrieved_max+1)
    precision_p = [0]*(retrieved_max+1)

    total_retrieved = len(retrieved_junctions)

    all_tp = 0
    all_t = 0
    all_p = 0

    # Count true positives and all true
    for j in true_junctions:
        recall_t[j[3]-1] += 1
        all_t += 1

        # If j is in retrieved_junctions, bisect_left returns either the correct
        #   index or index+1, depending on whether the coverage value is greater
        #   or less than the true value
        index = bisect_left(retrieved_junctions,j)
        if index == total_retrieved:
            continue

        # Check index
        match = retrieved_junctions[index]
        if j[2] == match[2] and j[1] == match[1] and j[0] == match[0]:
            recall_tp[j[3]-1] += 1
            precision_tp[match[3]-1] += 1
            all_tp += 1
        else:
            # Check index-1
            index -= 1

            if index < 0:
                continue

            match = retrieved_junctions[index]
            if j[2] == match[2] and j[1] == match[1] and j[0] == match[0]:
                recall_tp[j[3]-1] += 1
                precision_tp[match[3]-1] += 1
                all_tp += 1


    # Count all positives
    for j in retrieved_junctions:
        precision_p[j[3]-1] += 1
        all_p += 1


    print 'Overall recall = %f' % (float(all_tp)/all_t)
    print 'Overall precision = %f' % (float(all_tp)/all_p)

    # Construct CDFs of precision and recall for plots
    recall_cdf = [0] * len(recall_tp)


    #print recall_tp
    #print recall_t

    # keep track of current value of tp and p CDFs
    curr_recall_tp = 0
    curr_recall_t = 0
    for i in xrange(len(recall_tp)):
        curr_recall_tp += recall_tp[i]
        curr_recall_t += recall_t[i]
        recall_cdf[i] = float(curr_recall_tp) / curr_recall_t

    precision_cdf = [0] * len(precision_tp)
    # keep track of current value of tp and t CDFs
    curr_precision_tp = 0
    curr_precision_p = 0
    for i in xrange(len(precision_tp)):
        curr_precision_tp += precision_tp[i]
        curr_precision_p += precision_p[i]
        precision_cdf[i] = float(curr_precision_tp) / curr_precision_p

    print recall_cdf
    print precision_cdf

    # Graph
    plt.figure()
    plt.plot(range(1,len(recall_cdf)+1), recall_cdf)
    plt.xlabel('Coverage')
    plt.ylabel('Recall')
    plt.savefig('recall.png')

    plt.figure()
    plt.plot(range(1,len(precision_cdf)+1), precision_cdf)
    plt.xlabel('Coverage')
    plt.ylabel('Precision')
    plt.savefig('precision.png')


if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--true-junctions-bed', type=str, required=True, 
        help='Full path of BED file containing true junctions')
    parser.add_argument('--retrieved-junctions-bed', type=str, required=True, 
        help='Full path of BED file containing junctions retrieved by aligner')

    args = parser.parse_args(sys.argv[1:])
    
    with open(args.true_junctions_bed) as true_junctions_bed_stream:
        true_junctions = true_junctions_from_bed(true_junctions_bed_stream)
    
    with open(args.retrieved_junctions_bed) as retrieved_junctions_bed_stream:
        retrieved_junctions = retrieved_junctions_from_bed(retrieved_junctions_bed_stream)

    compare_junctions(true_junctions, retrieved_junctions)
