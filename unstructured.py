#!/usr/bin/env python
# coding: utf-8


from sarscov2_util import *



aln_file_1 = 'alignments/RR_nCov_alignment2_021120_muscle.fa'


# Gets secondary structure IN THE REVERSE COMPLEMENT
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
    betacov_intervals = get_ref_intervals_from_file(aln_file_1, cutoff=0.9, MIN_SIZE=15)
    unpaired_intervals = get_unpaired_intervals(per_position_unpaired, unpaired_cutoff=0.5, min_size=15)
    overlap_intervals = get_interval_overlap_size(unpaired_intervals, betacov_intervals, min_size=15)

    overlap_intervals = np.array(overlap_intervals)
    p_vals = np.array([x[3] for x in overlap_intervals])
    top_overlaps = overlap_intervals[np.argsort(p_vals)[::-1]]
    csv_file = 'results/unpaired_betacov_overlaps.csv'
    f = open(csv_file, 'w')
    for x in top_overlaps:
        f.write('%d,%d,%f,%f,%s\n' % (int(x[0]), int(x[1]-1), float(x[2]), float(x[3]), ref_seq[(int(x[0])-1):int(x[1]-1)]))
    f.close()
    print(len(unpaired_intervals))
    print(len(overlap_intervals))


def write_unpaired_sarscov2(per_position_unpaired, ref_seq):
    sarscov2_intervals = get_ref_intervals_from_file("alignments/gisaid_mafft_ncbi.fa", cutoff=0.97, MIN_SIZE=15)
    unpaired_intervals = get_unpaired_intervals(per_position_unpaired, unpaired_cutoff=0.6, min_size=15)
    overlap_intervals = get_interval_overlap_size(unpaired_intervals, sarscov2_intervals, min_size=15)

    overlap_intervals = np.array(overlap_intervals)
    p_vals = np.array([x[3] for x in overlap_intervals])
    top_overlaps = overlap_intervals[np.argsort(p_vals)[::-1]]
    csv_file = 'results/sarscov2_conserved_unstructured.csv'
    f = open(csv_file, 'w')
    for ii, x in enumerate(top_overlaps):
        f.write('SARS-CoV-2-conserved-unstructured-%d,%d-%d,%f,%f,%s\n' % (ii+1,int(x[0]), int(x[1]-1), float(x[2]), float(x[3]), ref_seq[(int(x[0])-1):int(x[1]-1)]))
    f.close()
    print(len(unpaired_intervals))
    print(len(overlap_intervals))


def write_unpaired_intervals(per_position_unpaired):
    unpaired_intervals = get_unpaired_intervals(per_position_unpaired, unpaired_cutoff=0.6, min_size=15)
    print(len(unpaired_intervals))
    # Print intervals in rank order of MCC secondary structure predictions
    f = open('results/unpaired_intervals_0.6_15.bed', 'w')
    for interval in unpaired_intervals:
        f.write('NC_045512.2\t%d\t%d\tint\t0\t.\n' % (interval[0], interval[1]))
    f.close()

if __name__ == "__main__":
    print("Getting and recording per-position unpaired")
    sequences = get_sequences(aln_file_1)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)
    per_position_unpaired = get_per_position_unpaired(ref_seq)


    f = open("results/probs_unpaired.txt", 'w')
    for cur_val in per_position_unpaired:
        f.write("%f\n" % cur_val)
    f.close()


    print("Recording unpaired intervals conserved in SARS-related sequences")
    write_unpaired_betacov(per_position_unpaired, ref_seq)
    print("Recording unpaired intervals conserved in SARS-CoV-2 sequences")
    write_unpaired_sarscov2(per_position_unpaired, ref_seq)
    print("Recording unpaired intervals")
    write_unpaired_intervals(per_position_unpaired)
