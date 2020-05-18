#!/usr/bin/env python
# coding: utf-8

from sarscov2_util import *
from matplotlib import pyplot as plt
from scipy import stats
import argparse


parser = argparse.ArgumentParser(description='Get conserved intervals, MEA structures around conserved intervals, and run relevant significance tests. Output file results/sars_related_conserved.csv has intervals conserved in SARS-related sequences and SARS-CoV-2 sequences to the desired conservation level. Output files results/conserved_intervals... have SARS-related conserved interval sequences, and results/intervals...bed has SARS-related conserved interval positions.')
parser.add_argument('--sarsr_aln_file', default='alignments/RR_nCov_alignment2_021120_muscle.fa', help="Alignment file for SARS-related sequences in fasta format")
parser.add_argument('--sarsr_aln_tag', default='RR_nCov_alignment2_021120_muscle', help="Tag for labeling results files for SARS-related conserved intervals")
parser.add_argument('--sarscov2_aln_file', default='alignments/gisaid_mafft_ncbi.fa', help="Alignment file for SARS-CoV-2 sequences in fasta format")
parser.add_argument('--min_len', type=int, default=15, help="Minimum length of reported conserved intervals")
parser.add_argument('--sarsr_conservation', type=float, default=1, help="Required conservation level for intervals in SARS-related sequences (float from 0 to 1)")
parser.add_argument('--sarscov2_conservation', type=float, default=0.99, help="Required conservation level for intervals in SARS-CoV-2 sequences (float from 0 to 1)")
parser.add_argument('--do_significance_tests', dest='do_significance_tests', action='store_true', help="If specified, will do significance tests for the overlap of SARS-related conserved intervals and SARS-CoV-2 conserved intervals")
parser.set_defaults(do_significance_tests=False)
args = parser.parse_args()

# ## SARS-related alignment analysis

# Alignment files: 
aln_file_1 = args.sarsr_aln_file
aln_file_2 = 'alignments/blast_alignment.fa'
aln_file_3 = 'alignments/BetaCov_comp_c990_full_nopA_clustalo_021920.fa'
aln_tag_1 = args.sarsr_aln_tag

sarscov2_aln_file1 = args.sarscov2_aln_file
sarscov2_aln_file2 = "alignments/gisaid_mafft_ncbi.fa"

MIN_SIZE = args.min_len
SARSR_CONSERVATION = args.sarsr_conservation
SARSCOV2_CONSERVATION = args.sarscov2_conservation
do_significance_tests = args.do_significance_tests

def print_conservation_datasets(aln_vals, ref_seq, all_intervals, aln_tag):
    # Print out the reference sequence
    f = open('results/refseq.txt', 'w')
    f.write("%s\n" % ref_seq)
    f.close()

    # Print intervals in rank order of MCC secondary structure predictions
    f = open('results/intervals_full_' + aln_tag + '.csv', 'w')
    for interval in all_intervals:
        f.write('%d,%d,%s,%s,%s,%f\n' % (interval.int_start - 19,interval.int_end + 20, interval.full_seq.replace('T', 'U'),                                      interval.secstruct, interval.seq, interval.mcc))
    f.close()

    # Print intervals in rank order of MCC secondary structure predictions, various formats
    f = open('results/intervals_full_' + aln_tag + '.csv', 'w')
    for interval in all_intervals:
        f.write('%d,%d,%s,%s,%s,%f\n' % (interval.int_start - 19,interval.int_end + 20, interval.full_seq.replace('T', 'U'),                                      interval.secstruct, interval.seq, interval.mcc))
    f.close()

    f = open('results/intervals_full_' + aln_tag + '.txt', 'w')
    for interval in all_intervals:
        f.write('%s %d %d\n' % (interval.seq, interval.int_start, interval.int_end-1))
        f.write('%s %f\n' % (interval.secstruct, interval.mcc))
        f.write('%s\n' % interval.full_seq.replace('T', 'U'))
        f.write('%s\n' % interval.conserved_region)
    f.close()

    f = open('results/intervals_full_' + aln_tag + '.bed', 'w')
    for interval in all_intervals:
        f.write('NC_045512.2\t%d\t%d\tint\t0\t.\n' % (interval.int_start, interval.int_end-1))
    f.close()

    f = open('results/conserved_intervals_' + aln_tag + '.fa', 'w')
    for interval in all_intervals:
        f.write('%s\n' % (ref_seq[interval.int_start:interval.int_end + 1]))
    f.close()

    # Print percentage conservation in alignment
    f = open('results/perc_conserved_' + aln_tag + '.txt', 'w')
    for aln_val in aln_vals:
        f.write('%f\n' % aln_val)
    f.close()



def get_structured_region_overlap_pval(aln_filename, cutoff, min_size):
    intervals = get_ref_intervals_from_file(aln_filename, cutoff=cutoff, MIN_SIZE=min_size)
    for region_key in regions.keys():
        region = regions[region_key]
        if get_interval_overlap_single(region, intervals):
            print(region)
    overlaps_list = get_perc_overlaps_rnd_trials([regions[region_key] for region_key in regions.keys()],                                                  intervals, len(ref_seq))
    overlaps_list = np.array(overlaps_list)
    print("P-value for overlap of structured regions: %f" % np.prod(overlaps_list))


# ## SARS-CoV-2 alignment analysis

def get_sarsr_sarscov2_conserved_intervals(sarsr_aln_file, sarscov2_aln_file, sarsr_cutoff, sarscov2_cutoff, min_size):
    sequences = get_sequences(sarscov2_aln_file)
    aln_counts = get_sequence_logo(sequences)
    (full_ref_seq, refseq) = get_ref_seq(sequences)
    aln_vals = get_aln_percs(aln_counts, full_ref_seq)

    intervals = get_ref_intervals_from_file(sarsr_aln_file, cutoff=sarsr_cutoff, MIN_SIZE=min_size)
    min_cons = []
    mean_cons = []
    interval_data = []
    for interval in intervals:
        cur_min_cons = min(aln_vals[interval[0]:interval[1]])
        min_cons += [cur_min_cons]
        mean_cons += [np.mean(np.array(aln_vals[interval[0]:interval[1]]))]
        if cur_min_cons > sarscov2_cutoff:
            interval_data += [(interval[0], interval[1], cur_min_cons)]
    return (interval_data, intervals, min_cons, mean_cons, aln_vals, refseq)

def sarsr_sarscov2_rnd_trials(sarsr_aln_file, sarscov2_aln_file, sarsr_cutoff, sarscov2_cutoff, min_size):
    (_, intervals, min_cons, mean_cons, aln_vals, refseq) = \
        get_sarsr_sarscov2_conserved_intervals(sarsr_aln_file, sarscov2_aln_file, \
            sarsr_cutoff, sarscov2_cutoff, min_size)

    rnd_min_cons = []
    rnd_mean_cons = []
    num_cons_wins = 0
    for ii in range(10000):
        rnd_min_cons = []
        rnd_mean_cons = []
        for interval in intervals:
            length = interval[1] - interval[0]
            start_idx = randint(0, len(refseq) - length)
            end_idx = start_idx + length
            rnd_min_cons += [min(aln_vals[start_idx:end_idx])]
            rnd_mean_cons += [np.mean(np.array(aln_vals[start_idx:end_idx]))]
        if np.sum(np.array(mean_cons) > sarscov2_cutoff) > np.sum(np.array(rnd_mean_cons) > sarscov2_cutoff):
            num_cons_wins += 1
    print("Number of trials with more SARSr conserved intervals 99%% conserved than random intervals: %d/10000" % num_cons_wins)
    print("Binomial test p-value for this difference: %.2E" % stats.binom_test(num_cons_wins, 10000))
    plt.hist([mean_cons, rnd_mean_cons], range=(0.95, 1), bins=10)
    plt.legend(["Conserved intervals", "Random intervals"])
    plt.show()

if __name__ == '__main__':
    print("Getting SARSr conservation intervals, sequence conservation percentages, and MEA structure predictions")
    sequences = get_sequences(aln_file_1)
    bcov_aln_counts = get_sequence_logo(sequences)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)
    aln_vals = get_aln_percs(bcov_aln_counts, full_ref_seq)
    ref_intervals = get_intervals_refseq(bcov_aln_counts, full_ref_seq, sequences, \
        cutoff=SARSR_CONSERVATION, MIN_SIZE=MIN_SIZE)


    all_intervals = get_all_secstructs_mea(ref_intervals, ref_seq, secstruct_interval=20)
    print("Printing SARSr conservation intervals, sequence conservation percentages, and MEA structure predictions")
    print_conservation_datasets(aln_vals, ref_seq, all_intervals, aln_tag_1)

    if do_significance_tests:
        print("Structured regions that overlap with SARSr-conserved regions")
        get_structured_region_overlap_pval(aln_file_1, SARSR_CONSERVATION, MIN_SIZE)

        print("Comparing SARSr sequence conservation in SARS-CoV-2 to random sequences using alignment: ")
        print("NCBI")
        sarsr_sarscov2_rnd_trials(aln_file_1, sarscov2_aln_file1, \
            SARSR_CONSERVATION, SARSCOV2_CONSERVATION, MIN_SIZE)
        print("GSAID")
        sarsr_sarscov2_rnd_trials(aln_file_1, sarscov2_aln_file2, \
            SARSR_CONSERVATION, SARSCOV2_CONSERVATION, MIN_SIZE)


    # Record SARS-related, SARS-CoV-2 conserved sequences
    print("Recording SARS-related, SARS-CoV-2 conserved sequences")
    (interval_data, intervals, min_cons, mean_cons, aln_vals, refseq) = \
        get_sarsr_sarscov2_conserved_intervals(aln_file_1, sarscov2_aln_file2, \
            SARSR_CONSERVATION, SARSCOV2_CONSERVATION, MIN_SIZE)

    interval_data = np.array(interval_data)
    idxs = np.argsort(np.array([x[2] for x in interval_data]))
    interval_data = interval_data[idxs[::-1]]
    f = open('results/sars_related_conserved.csv','w')
    for ii, interval in enumerate(interval_data):
        f.write("SARS-related-conserved-%d,%d-%d,%s,%f\n" % (ii+1,interval[0], interval[1], refseq[int(interval[0]-1):int(interval[1])],interval[2]))
    f.close()
