#!/usr/bin/env python
# coding: utf-8


from sarscov2_util import *
import os
import argparse

parser = argparse.ArgumentParser(description='Get conserved, unpaired stretches using Contrafold 2.0 base-pairing probabilities and RNAplfold. Output file results/unpaired_betacov_overlaps.csv has intervals that are unpaired at the specified thresholds and conserved in SARS-related sequences. Output file results/sarscov2_conserved_unstructured.csv has intervals that are unpaired and conserved in SARS-CoV-2 sequences. Output file results/unpaired_intervals.bed has all unpaired intervals at the specified threshold regardless of conservation. Output file results/rnaplfold_intervals.csv has unpaired intervals that overlap between Contrafold 2.0 and RNAplfold.')
parser.add_argument('--sarsr_aln_file', default='alignments/RR_nCov_alignment2_021120_muscle.fa', help="Alignment file for SARS-related sequences in fasta format")
parser.add_argument('--sarscov2_aln_file', default='alignments/gisaid_mafft_ncbi.fa', help="Alignment file for SARS-CoV-2 sequences in fasta format")
parser.add_argument('--min_len', type=int, default=15, help="Minimum length of reported conserved unpaired intervals")
parser.add_argument('--sarsr_conservation', type=float, default=0.9, help="Required conservation level for intervals in SARS-related sequences (float from 0 to 1)")
parser.add_argument('--sarscov2_conservation', type=float, default=0.97, help="Required conservation level for intervals in SARS-CoV-2 sequences (float from 0 to 1)")
parser.add_argument('--min_unpaired_prob', type=float, default=0.6, help="Minimum probability unpaired for the position to be classified as unstructured; all positions in unpaired intervals must have unpaired probabilities above this threshold.")
parser.add_argument('--do_significance_tests', dest='do_significance_tests', action='store_true', help="If specified, will do significance tests for the overlap of SARS-related conserved intervals and SARS-CoV-2 conserved intervals")

args = parser.parse_args()

aln_file_1 = args.sarsr_aln_file
sarscov2_aln_file = args.sarscov2_aln_file
MIN_SIZE = args.min_len
SARSR_CONSERVATION = args.sarsr_conservation
SARSCOV2_CONSERVATION = args.sarscov2_conservation
MIN_UNPAIRED = args.min_unpaired_prob
do_significance_tests = args.do_significance_tests


def get_unpaired_probs(int_start, int_end, ref_seq):
    sequence = ref_seq[int_start:int_end]
    bp_matrix = bpps(sequence, package='contrafold_2')
    paired_probs = np.sum(bp_matrix, axis=0)
    return 1 - paired_probs

# Just as in RNAz, use 120 bp windows and slide 40 bp
def get_per_position_unpaired(ref_seq, window_size=120, slide=40, print_freq=1000):
    total_unpaired_probs = np.zeros(len(ref_seq))
    num_unpaired_probs = np.zeros(len(ref_seq))
    ii = 0
    while ii < len(ref_seq):
        window_len = min(window_size, len(ref_seq) - ii)
        if ii % print_freq == 0:
            print("Window: %d to %d" % (ii, ii + window_len))
        unpaired_probs = get_unpaired_probs(ii, ii + window_len, ref_seq)
        total_unpaired_probs[ii:(ii + window_len)] += unpaired_probs
        num_unpaired_probs[ii:(ii + window_len)] += np.ones(window_len)
        ii += slide
    return total_unpaired_probs/num_unpaired_probs


def get_unpaired_intervals(probs_unpaired, unpaired_cutoff=0.9, min_size=15):
    intervals = []
    start_int = -1
    end_int = -1
    in_int = False
    int_size = 0
    for ii in range(len(probs_unpaired)):
        if (probs_unpaired[ii] < unpaired_cutoff):
            if in_int:
                if int_size > min_size:
                    avg_unpaired = np.mean(np.array(probs_unpaired[start_int:end_int]))
                    min_unpaired = min(probs_unpaired[start_int:end_int])
                    intervals += [(start_int + 1, end_int, avg_unpaired, min_unpaired)]
                int_size = 0
                in_int = False
        else:
            if not in_int:
                in_int = True
                start_int = ii
            int_size += 1
            end_int = ii
    if in_int:
        if int_size > min_size:
            avg_unpaired = np.mean(np.array(probs_unpaired[start_int:end_int]))
            min_unpaired = min(probs_unpaired[start_int:end_int])
            intervals += [(start_int + 1, end_int, avg_unpaired, min_unpaired)]
    return intervals


def write_unpaired_betacov(per_position_unpaired, ref_seq):
    betacov_intervals = get_ref_intervals_from_file(aln_file_1, cutoff=SARSR_CONSERVATION, MIN_SIZE=MIN_SIZE)
    unpaired_intervals = get_unpaired_intervals(per_position_unpaired, unpaired_cutoff=MIN_UNPAIRED, min_size=MIN_SIZE)
    overlap_intervals = get_interval_overlap_size(unpaired_intervals, betacov_intervals, min_size=MIN_SIZE)

    overlap_intervals = np.array(overlap_intervals)
    p_vals = np.array([x[3] for x in overlap_intervals])
    top_overlaps = overlap_intervals[np.argsort(p_vals)[::-1]]
    csv_file = 'results/unpaired_betacov_overlaps.csv'
    f = open(csv_file, 'w')
    f.write('%s,%s,%s,%s,%s\n' % ("Start index", "end index", "average unpaired probability in interval", "minimum unpaired probability in interval", "interval sequence"))
    for x in top_overlaps:
        f.write('%d,%d,%f,%f,%s\n' % (int(x[0]), int(x[1]-1), float(x[2]), float(x[3]), ref_seq[(int(x[0])-1):int(x[1]-1)]))
    f.close()


def write_unpaired_sarscov2(per_position_unpaired, ref_seq):
    sarscov2_intervals = get_ref_intervals_from_file(sarscov2_aln_file, cutoff=SARSCOV2_CONSERVATION, MIN_SIZE=MIN_SIZE)
    unpaired_intervals = get_unpaired_intervals(per_position_unpaired, unpaired_cutoff=MIN_UNPAIRED, min_size=MIN_SIZE)
    overlap_intervals = get_interval_overlap_size(unpaired_intervals, sarscov2_intervals, min_size=MIN_SIZE)

    overlap_intervals = np.array(overlap_intervals)
    p_vals = np.array([x[3] for x in overlap_intervals])
    top_overlaps = overlap_intervals[np.argsort(p_vals)[::-1]]
    csv_file = 'results/sarscov2_conserved_unstructured.csv'
    f = open(csv_file, 'w')
    f.write('%s,%s,%s,%s,%s\n' % ("Start index", "end index", "average unpaired probability in interval", "minimum unpaired probability in interval", "interval sequence"))
    for ii, x in enumerate(top_overlaps):
        f.write('SARS-CoV-2-conserved-unstructured-%d,%d-%d,%f,%f,%s\n' % (ii+1,int(x[0]), int(x[1]-1), float(x[2]), float(x[3]), ref_seq[(int(x[0])-1):int(x[1]-1)]))
    f.close()

def write_unpaired_intervals(unpaired_intervals):
    # Print intervals in rank order of MCC secondary structure predictions
    f = open('results/unpaired_intervals.bed', 'w')
    for interval in unpaired_intervals:
        f.write('NC_045512.2\t%d\t%d\tint\t0\t.\n' % (interval[0], interval[1]))
    f.close()

def write_rnaplfold_intervals(intervals):
    f = open('results/rnaplfold_intervals.csv', 'w')
    f.write('%s,%s\n' % ("Start index", "end index"))
    for interval in intervals:
        f.write("%d, %d\n" % (interval[0], interval[1]))
    f.close()

def get_rnaplfold_intervals(RNAplfold_filename, int_size=15, cutoff=0.25):
    f = open(RNAplfold_filename)
    rnaplfold_lines = f.readlines()
    f.close()

    rnaplfold_intervals = []
    for rnaplfold_line in rnaplfold_lines[(int_size+2):]:
        rnaplfold_items = rnaplfold_line.split('\t')
        int_start = int(rnaplfold_items[0])
        if float(rnaplfold_items[-1]) > cutoff:
            last_int = (-1, -1)
            if len(rnaplfold_intervals) > 0:
                last_int = rnaplfold_intervals[-1]
            # Continuation of last interval or not?
            if int_start < last_int[1]:
                rnaplfold_intervals = rnaplfold_intervals[:-1]
                rnaplfold_intervals += [(last_int[0], int_start + int_size)]
            else:
                rnaplfold_intervals += [(int_start, int_start+int_size)]
    print([x[0] for x in rnaplfold_intervals])

    return rnaplfold_intervals

def get_rnaplfold_per_position_unpaired(RNAplfold_filename):
    f = open(RNAplfold_filename)
    rnaplfold_lines = f.readlines()
    f.close()

    rnaplfold_unpaired_probabilities = []
    for rnaplfold_line in rnaplfold_lines[2:]:
        rnaplfold_unpaired_probabilities += [float(rnaplfold_line.split('\t')[1])]

    return rnaplfold_unpaired_probabilities

if __name__ == "__main__":
    print("Getting and recording per-position unpaired")
    sequences = get_sequences(aln_file_1)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)

    probs_unpaired_file = "results/probs_unpaired.txt"
    per_position_unpaired = []
    if os.path.isfile(probs_unpaired_file):
        f = open(probs_unpaired_file)
        per_position_unpaired = [float(x) for x in f.readlines()]
        f.close()
    else: 
        per_position_unpaired = get_per_position_unpaired(ref_seq)
        f = open(probs_unpaired_file, 'w')
        for cur_val in per_position_unpaired:
            f.write("%f\n" % cur_val)
        f.close()

    unpaired_intervals = get_unpaired_intervals(per_position_unpaired, unpaired_cutoff=MIN_UNPAIRED, min_size=MIN_SIZE)
    rnaplfold_unpaired_probabilities = get_rnaplfold_per_position_unpaired("RNAplfold/SARSCOV2_lunp")
    rnaplfold_intervals = get_unpaired_intervals(rnaplfold_unpaired_probabilities, unpaired_cutoff=MIN_UNPAIRED, min_size=MIN_SIZE)
    overlap_intervals = get_interval_overlap_size(unpaired_intervals, rnaplfold_intervals, min_size=5)
    print("%d/%d Contrafold 2.0 intervals overlap with RNAplfold" % (len(overlap_intervals), len(unpaired_intervals)))
    if do_significance_tests:
        num_overlaps = get_num_overlaps_rnd_trials_size(unpaired_intervals, rnaplfold_intervals, len(ref_seq), min_size=5)
        print("P-value for overlap: %.2E\n" % (np.sum(np.array(num_overlaps) >= len(overlap_intervals))/len(num_overlaps)))
    overlap_intervals = get_interval_overlap_size(rnaplfold_intervals, unpaired_intervals, min_size=5)
    print("%d/%d RNAplfold intervals overlap with Contrafold 2.0" % (len(overlap_intervals), len(rnaplfold_intervals)))
    if do_significance_tests:
        num_overlaps = get_num_overlaps_rnd_trials_size(rnaplfold_intervals, unpaired_intervals, len(ref_seq), min_size=5)
        print("P-value for overlap: %.2E\n" % (np.sum(np.array(num_overlaps) >= len(overlap_intervals))/len(num_overlaps)))
    
    print("Recording RNAplfold intervals that overlap with Contrafold 2.0 intervals")
    write_rnaplfold_intervals(overlap_intervals)
    print("Recording unpaired intervals conserved in SARS-related sequences")
    write_unpaired_betacov(per_position_unpaired, ref_seq)
    print("Recording unpaired intervals conserved in SARS-CoV-2 sequences")
    write_unpaired_sarscov2(per_position_unpaired, ref_seq)
    print("Recording unpaired intervals")
    write_unpaired_intervals(unpaired_intervals)
