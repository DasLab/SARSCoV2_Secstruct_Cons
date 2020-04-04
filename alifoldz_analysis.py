#!/usr/bin/env python
# coding: utf-8


from sarscov2_util import *
import os
from Bio import AlignIO
from scipy import stats



def make_stockholm_alignment(window_start, window_end, sequences, alignment_file):
    # Make fasta file alignment with desired range
    f = open('alignments/tmp_alignment.fa', 'w')
    
    for sequence in sequences:
        f.write('> %s\n' % sequence.name)
        cur_seq = sequence.seq[window_start:window_end]
        ii = 0
        while ii < len(cur_seq):
            f.write("%s\n" % cur_seq[ii:(ii+60)])
            ii += 60
    f.close()

    # Make stockholm alignment file with desired range
    align = AlignIO.read('alignments/tmp_alignment.fa', 'fasta')
    AlignIO.write(align, alignment_file, 'stockholm')

    os.remove('alignments/tmp_alignment.fa')



def make_fa_alignment(window_start, window_end, sequences, alignment_file):
    # Make fasta file alignment with desired range
    f = open(alignment_file, 'w')
    
    for sequence in sequences:
        f.write('> %s\n' % sequence.name)
        cur_seq = sequence.seq[window_start:window_end]
        ii = 0
        while ii < len(cur_seq):
            f.write("%s\n" % cur_seq[ii:(ii+60)])
            ii += 60
    f.close()


if __name__=='__main__':
    sequences = get_sequences('alignments/betacorona-genome-ref-oc43.fa')
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)
    
    # Make windowed alignments for use with rscape
    ii = 0
    window_id = 0
    while ii < len(full_ref_seq):
        window_start = ii
        window_end = ii + 120
        alignment_file = 'rscape/windows_betacov/alignment_' + str(window_id) + '.sto'
        make_stockholm_alignment(window_start, window_end, sequences, alignment_file)
        ii += 40
        window_id += 1
    
    sequences = get_sequences('alignments/RR_nCov_alignment2_021120_muscle.fa')
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)

    # Make windowed alignments for use with alifoldz
    ii = 0
    window_id = 0
    window_id_dict = {}
    while ii < len(full_ref_seq):
        window_start = ii
        window_end = ii + 120
        alignment_file = 'alifoldz/windows/alignment_' + str(window_id) + '.fa'
        make_fa_alignment(window_start, window_end, sequences, alignment_file)
        window_id_dict[window_id] = (window_start, window_end)
        ii += 40
        window_id += 1
    
    
    
    # Determine alifoldz cutoff
    # Compare alifoldz windows to rnaz windows:
    f = open("alifoldz/alifoldz_shuffle_data.dat")
    alifoldz_lines = f.readlines()
    f.close()
    
    # Chosen arbitrarily for now, rechoose based on shuffled sequences.
    alifoldz_vals = []
    for alifoldz_line in alifoldz_lines:
        alifoldz_items = alifoldz_line.strip('\n').split(' ')
        if alifoldz_items[3] == '':
            continue
        alifoldz_vals += [float(alifoldz_items[3])]
    alifoldz_vals = np.array(alifoldz_vals)
    CUTOFF = np.quantile(alifoldz_vals, 0.01)
    print(CUTOFF)
    
    
    
    aln_file_1 = 'alignments/RR_nCov_alignment2_021120_muscle.fa'
    rnaz_windows_f = 'rnaz_data/rnaz_RR_nCov_alignment2_021120_muscle.out'
    rnaz_windows_dict = get_rnaz_windows(rnaz_windows_f, aln_file_1)
    
    
    
    # Compare alifoldz windows to rnaz windows:
    f = open("alifoldz/alifoldz_data.dat")
    alifoldz_lines = f.readlines()
    f.close()
    
    # Do alifoldz intervals with z scores below cutoff based on shuffled alifoldz
    # overlap with RNAz intervals?
    alifoldz_intervals = []
    rnaz_overlaps = 0
    num_cutoff = 0
    f = open("results/rnaz_alifoldz_windows.csv", 'w')
    for ii, alifoldz_line in enumerate(alifoldz_lines):
        alifoldz_items = alifoldz_line.strip('\n').split(' ')
        if alifoldz_items[3] == '':
            continue
        if float(alifoldz_items[3]) <= CUTOFF:
            num_cutoff += 1
            (int_start, int_end) = window_id_dict[ii]
            ref_int_start = get_position(full_ref_seq, int_start)
            ref_int_end = get_position(full_ref_seq, int_end) - 1
            alifoldz_intervals += [(ref_int_start, ref_int_end)]
            if (ref_int_start, ref_int_end) in rnaz_windows_dict.keys():
                rnaz_overlaps += 1
                (seq, rnaz_secstruct, rnaz_z, rnaz_p) = rnaz_windows_dict[(ref_int_start, ref_int_end)]
                alifoldz_z = float(alifoldz_items[3])
                f.write('%d,%d,%s,%s,%.2f,%.4f,%.2f\n' % (ref_int_start, ref_int_end, seq, rnaz_secstruct, rnaz_z, rnaz_p, alifoldz_z))
    f.close()
    
    total_rnaz = len(rnaz_windows_dict.keys())
    total_alifoldz = num_cutoff
    total_windows = len(window_id_dict.keys())
    total_overlap = rnaz_overlaps
    print("Total overlapping intervals: %d; Total RNAz intervals: %d; Total alifoldz intervals: %d; Total Windows: %d\n" %      (total_overlap, total_rnaz, total_alifoldz, total_windows))
    rv = stats.hypergeom(total_windows, total_alifoldz, total_rnaz)
    p_value = 1 - rv.cdf(total_overlap)
    print("P value for overlap: %.2E" % p_value)
    
    # Do alifoldz intervals overlap with structured elements?
    print("Structured elements that overlap with alifoldz windows:")
    for region_key in regions.keys():
        region_interval = regions[region_key]
        region_start = get_position_full(full_ref_seq, region_interval[0])
        region_end = get_position_full(full_ref_seq,  region_interval[1])
        if get_interval_overlap_single([region_start, region_end], alifoldz_intervals):
            print(region_key)