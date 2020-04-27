#!/usr/bin/env python
# coding: utf-8


from sarscov2_util import *
from matplotlib import pyplot as plt



aln_file_1 = 'alignments/RR_nCov_alignment2_021120_muscle.fa'
rnaz_results_1 = 'rnaz_data/results_sorted_RR_nCov_alignment2_021120_muscle.dat'
scanfold_file = 'scanfold_data/tables1.csv'

def get_rnaz_intervals(rnaz_filename):
    f = open(rnaz_filename)
    rnaz_lines = f.readlines()
    f.close()
    
    intervals = []
    
    for rnaz_line in rnaz_lines:
        print(rnaz_line[0:6])
        if len(rnaz_line) > 5 and rnaz_line[0:5] == 'locus':
            rnaz_items = rnaz_line.split('\t')
            start_idx = int(rnaz_items[2])
            end_idx = int(rnaz_items[3])
            intervals += [(start_idx, end_idx)]
            
    return intervals


def get_rnaz_intervals_windows(rnaz_filename):
    f = open(rnaz_filename)
    rnaz_lines = f.readlines()
    f.close()
    
    intervals = []
    
    for rnaz_line in rnaz_lines:
        if len(rnaz_line) > 6 and rnaz_line[0:6] == 'window':
            rnaz_items = rnaz_line.split('\t')
            start_idx = int(rnaz_items[3])
            end_idx = int(rnaz_items[4])
            p_val = float(rnaz_items[-1])
            intervals += [(start_idx, end_idx, p_val)]
            
    return intervals


def convert_rnaz_intervals(rnaz_intervals, aln_file): 
    sequences = get_sequences(aln_file)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)
    return convert_intervals_refseq(rnaz_intervals, full_ref_seq)

def get_rnaz_intervals_refseq(rnaz_filename, aln_filename):
    rnaz_intervals = get_rnaz_intervals_windows(rnaz_filename)
    return convert_rnaz_intervals(rnaz_intervals, aln_filename)


# Write all RNAz windows and loci
def write_rnaz_windows_loci():
    f = open(rnaz_results_1)
    rnaz_lines = f.readlines()
    f.close()

    sequences = get_sequences(aln_file_1)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)

    f = open('results/rnaz_loci.csv', 'w')
    for rnaz_line in rnaz_lines:
        if len(rnaz_line) > 6 and rnaz_line[0:6] == 'window':
            rnaz_items = rnaz_line.split('\t')
            start_idx = int(rnaz_items[3])
            start_idx = get_position(full_ref_seq, start_idx)
            end_idx = int(rnaz_items[4])
            end_idx = get_position(full_ref_seq, end_idx)
            z_score = float(rnaz_items[-4])
            sci = float(rnaz_items[-3])
            p_val = float(rnaz_items[-1])
            f.write('%s,%d-%d,%f,%f,%f\n' % (rnaz_items[0], start_idx, end_idx, z_score, sci, p_val))
        if len(rnaz_line) > 5 and rnaz_line[0:5] == 'locus':
            rnaz_items = rnaz_line.split('\t')
            start_idx = int(rnaz_items[2])
            start_idx = get_position(full_ref_seq, start_idx)
            end_idx = int(rnaz_items[3])
            end_idx = get_position(full_ref_seq, end_idx)
            p_val = float(rnaz_items[-2])
            z = float(rnaz_items[-1])
            f.write("%s,%d-%d,%f,%f\n" % (rnaz_items[0], start_idx, end_idx, z, p_val))
    f.close()

# Write all RNAz windows, secondary structures, z scores to a file
# Save all RNAz windows to a dictionary
def write_rnaz_windows():
    rnaz_windows_f = 'rnaz_data/rnaz_RR_nCov_alignment2_021120_muscle.out'
    f = open(rnaz_windows_f)
    rnaz_windows = f.readlines()
    f.close()

    rnaz_windows_dict = {}

    sequences = get_sequences(aln_file_1)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)

    csv_file = 'results/rnaz_RR_nCov_alignment2_021120_muscle.csv'
    ii = 0
    f = open(csv_file, 'w')
    while ii < len(rnaz_windows):
        line = rnaz_windows[ii]
        if "RNA-class probability" in line:
            p_val = float(line.split(' ')[-1])
        if 'NC_045512.2' in line:
            int_start = int(line.split('/')[1].split('-')[0])
            int_start = get_position(full_ref_seq, int_start)
            int_end = int(line.split('/')[1].split('-')[1])-1
            int_end = get_position(full_ref_seq, int_end)
            ii += 1
            line = rnaz_windows[ii]
            seq = line.strip('\n').replace('-', '')
            ii += 1
            line = rnaz_windows[ii]
            secstruct = line.split(' ')[0].replace('-', '')
            z_score_attempt = line.split(' ')[6].strip(',')
            if z_score_attempt == '': # this is when z-score is positive
                z_score_attempt = line.split(' ')[7].strip(',')
            z_score = float(z_score_attempt)
            f.write('%d,%d,%s,%s,%f,%f\n' % (int_start, int_end, seq, secstruct, z_score, p_val))
            rnaz_windows_dict[(int_start, int_end)] = (seq, secstruct, z_score, p_val)
        ii += 1
    f.close()
    
    return rnaz_windows_dict

def get_rnaz_scanfold_overlap(rnaz_windows_dict, ref_seq, p_val_cutoff=0.05):
    f = open(scanfold_file)
    scanfold_lines = f.readlines()[1:]
    f.close()

    scanfold_intervals = []
    for scanfold_line in scanfold_lines:
        scanfold_items = scanfold_line.replace('\ufeff', '').split(',')
        int_ends = (int(scanfold_items[0]), int(scanfold_items[1]))
        p_val = float(scanfold_items[5])
        if p_val < p_val_cutoff:
            scanfold_intervals += [int_ends]

    rnaz_intervals = []
    for rnaz_key in rnaz_windows_dict.keys():
        rnaz_window = rnaz_windows_dict[rnaz_key]
        if rnaz_window[3] > 0.9:
            rnaz_intervals += [rnaz_key]
    
    print("Number of scanfold intervals: %d" % len(scanfold_intervals))
    print("Number of rnaz windows: %d" % len(rnaz_intervals))

    overlap_intervals = get_interval_overlap_size(rnaz_intervals, scanfold_intervals, full_int=True)
    print("Number of overlapping intervals: %d/%d" % (len(overlap_intervals), len(rnaz_intervals)))
    num_overlaps = get_num_overlaps_rnd_trials_size(rnaz_intervals, scanfold_intervals, len(ref_seq), full_int=True)
    print("P-value for overlap: %f\n" % (np.sum(np.array(num_overlaps) >= len(overlap_intervals))/len(num_overlaps)))
    
    overlap_intervals = get_interval_overlap(scanfold_intervals, rnaz_intervals)
    print("Number of overlapping intervals: %d/%d" % (len(overlap_intervals), len(scanfold_intervals)))
    num_overlaps = get_num_overlaps_rnd_trials(scanfold_intervals, rnaz_intervals, len(ref_seq))
    print("P-value for overlap: %f\n" % (np.sum(np.array(num_overlaps) >= len(overlap_intervals))/len(num_overlaps)))


# Write all SARS-CoV-2 overlapping structured intervals to a file
def write_structured_conserved_regions(rnaz_windows_dict):
    sarscov2_intervals = get_ref_intervals_from_file("alignments/gisaid_mafft_ncbi.fa", cutoff=0.97, MIN_SIZE=15)
    rnaz_intervals = get_rnaz_intervals_refseq(rnaz_results_1, aln_file_1)

    overlap_intervals = get_interval_overlap_size(rnaz_intervals, sarscov2_intervals, min_size=15)
    overlap_intervals = np.array(overlap_intervals)
    p_vals = np.array([x[2] for x in overlap_intervals])
    top_overlaps = overlap_intervals[np.argsort(p_vals)[::-1]]
    csv_file = 'results/sarscov2_conserved_structured.csv'
    f = open(csv_file, 'w')
    ii = 0
    for x in top_overlaps:
        interval = (int(x[0])-1, int(x[1])-2)
        if interval in rnaz_windows_dict.keys():
            ii += 1
            rnaz_window = rnaz_windows_dict[interval]
            f.write("SARS-CoV-2-conserved-structured-%d,%d-%d,%s,%s,%f,%f\n" % (ii, interval[0], interval[1], rnaz_window[0], rnaz_window[1], float(rnaz_window[2]), float(rnaz_window[3])))
    f.close()


def write_structured_sarsr_conserved_regions(ref_seq):
    betacov_intervals = get_ref_intervals_from_file(aln_file_1, cutoff=0.9, MIN_SIZE=15)
    rnaz_intervals = get_rnaz_intervals_refseq(rnaz_results_1, aln_file_1)

    overlap_intervals = get_interval_overlap_size(rnaz_intervals, betacov_intervals, min_size=15)
    overlap_intervals = np.array(overlap_intervals)
    p_vals = np.array([x[2] for x in overlap_intervals])
    top_overlaps = overlap_intervals[np.argsort(p_vals)[::-1]]

    csv_file = 'results/rnaz_betacov_overlaps.csv'
    f = open(csv_file, 'w')
    for x in top_overlaps:
        f.write("%d,%d,%f,%s\n" % (int(x[0]), int(x[1]-1), float(x[2]), ref_seq[(int(x[0])-1):int(x[1]-1)]))
    f.close()



def test_rnaz_interval_overlaps_pvals(ref_seq):
    sarscov2_intervals = get_ref_intervals_from_file("alignments/gisaid_mafft_ncbi.fa", cutoff=0.97, MIN_SIZE=15)
    rnaz_intervals = get_rnaz_intervals_refseq(rnaz_results_1, aln_file_1)
    overlap_intervals = get_interval_overlap(rnaz_intervals, sarscov2_intervals)
    print("Number of structured regions: %d" % len(rnaz_intervals))
    print("Number of conserved regions: %d\n" % len(sarscov2_intervals))

    # Number of canonical regions overlapping with RNAz structure predictions, and how unlikely is this?
    print("Canonical regions overlapping with RNAz structure predictions")
    overlap_intervals = get_interval_overlap([regions[region_key] for region_key in regions.keys()], rnaz_intervals)
    print("Number of overlapping intervals: %d/%d" % (len(overlap_intervals), len(regions.keys())))
    num_overlaps = get_num_overlaps_rnd_trials([regions[region_key] for region_key in regions.keys()], rnaz_intervals, len(ref_seq))
    print("P-value for overlap: %f\n" % (np.sum(np.array(num_overlaps) >= len(overlap_intervals))/len(num_overlaps)))

    # Number of  RNAz structure predictions overlapping with conserved intervals, and how unlikely is this? 
    print("RNAz structure predictions overlapping with conserved intervals")
    overlap_intervals = get_interval_overlap_size(rnaz_intervals, sarscov2_intervals, min_size=15)
    overlap_intervals = np.array(overlap_intervals)
    p_vals = np.array([x[2] for x in overlap_intervals])
    print("Number of overlapping intervals: %d/%d" % (len(overlap_intervals), len(rnaz_intervals)))
    num_overlaps = get_num_overlaps_rnd_trials_size(rnaz_intervals, sarscov2_intervals, len(ref_seq), min_size=15)
    print("P-value for overlap: %f\n" % (np.sum(np.array(num_overlaps) >= len(overlap_intervals))/len(num_overlaps)))

    # Number of conserved intervals overlapping with RNAz structure predictions, and how unlikely is this? 
    print("Conserved intervals overlapping with RNAz structure predictions")
    overlap_intervals = get_interval_overlap_size(sarscov2_intervals, rnaz_intervals, min_size=15)
    print("Number of overlapping intervals: %d/%d" % (len(overlap_intervals), len(sarscov2_intervals)))
    num_overlaps = get_num_overlaps_rnd_trials_size(sarscov2_intervals, rnaz_intervals, len(ref_seq), min_size=15)
    print("P-value for overlap: %f\n" % (np.sum(np.array(num_overlaps) >= len(overlap_intervals))/len(num_overlaps)))


# What percent of the genome is in a structured region?
def get_percent_structured(ref_seq):
    rnaz_intervals = get_rnaz_intervals_refseq(rnaz_results_1, aln_file_1)

    total_structured = 0
    for rnaz_interval in rnaz_intervals:
        total_structured += rnaz_interval[1] - rnaz_interval[0]
    return total_structured/len(ref_seq)

if __name__ == '__main__':
    sequences = get_sequences(aln_file_1)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)

    print("Writing RNAz loci and windows to .csv files\n")
    write_rnaz_windows_loci()
    rnaz_windows_dict = write_rnaz_windows()
    
    #print("Recording structured regions conserved in SARS-CoV-2\n")
    #write_structured_conserved_regions(rnaz_windows_dict)
    
    print("Getting overlaps between RNAz intervals and ScanFold data\n")
    get_rnaz_scanfold_overlap(rnaz_windows_dict, ref_seq)

    #print("Recording structured regions conserved in SARS-CoV-2\n")
   #
    #print("Recording structured regions conserved in SARS-related viruses\n")
    #write_structured_sarsr_conserved_regions(ref_seq)
#
    #print("Testing overlaps between RNAz intervals, known structures, and conserved intervals")
    #test_rnaz_interval_overlaps_pvals(ref_seq)
    #
    #perc_structured = get_percent_structured(ref_seq)
    #print("Percent of genome in an RNAz structured window: %f\n" % perc_structured)
