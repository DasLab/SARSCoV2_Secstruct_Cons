import numpy as np
from arnie.pfunc import pfunc
from arnie.bpps import bpps
import arnie.utils as utils
from arnie.mea.mea import MEA
from random import randint

class Seq:
    def __init__(self, seq='', name='', is_ncov=False, is_reference=False, is_shortlist=True):
        self.seq = seq
        self.name = name
        self.is_ncov = is_ncov
        self.is_reference = is_reference
        self.is_shortlist = is_shortlist
        self.is_new_shortlist = is_shortlist

        
class Interval:
    def __init__(self, seq='', int_start=-1, int_end=-1, secstruct='', mcc=-1, conserved_region='', full_seq=''):
        self.seq = seq
        self.int_start=int_start
        self.int_end=int_end
        self.secstruct=secstruct
        self.mcc=mcc
        self.conserved_region=conserved_region
        self.full_seq=full_seq

refseq_names = ['NC_045512', 'NC_045512.2', 'Wuhan-Hu-1']
# Positions of known Rfam and structured elements
regions = {'fse': (13475, 13542), 'sl23': (45, 75), 'sl5': (150, 295), 'pk': (29608, 29657), 's2m': (29723, 29773)} 
pos_idx = {'A': 0, 'G': 1, 'C': 2, 'T': 3, '-': 4}

# Collect sequences from alignment
def get_sequences(alignment_file):
    f = open(alignment_file)
    new_alignment_lines = f.readlines()
    f.close()

    sequences = []
    
    cur_tag = ''
    cur_seq = ''
    for cur_line in new_alignment_lines:
        if len(cur_line) > 0 and cur_line[0] == '>':
            if cur_tag != '':
                is_reference = False
                for refseq_name in refseq_names:
                    if refseq_name in cur_tag:
                        is_reference = True
                if ('Query' not in cur_tag) and \
                    (len(cur_seq.strip('-')) > 1000):
                    sequences += [Seq(seq=cur_seq, name=cur_tag, is_reference=is_reference)]
                cur_seq = ''
            # The following tag processing handles a variety of alignment tag formats
            cur_tag = cur_line.strip('\n').strip('>').split(' ')[0].split('|')[-1] 
        else:
            cur_seq += cur_line.strip('\n').upper()
    if ('Query' not in cur_tag) and \
        (len(cur_seq.strip('-')) > 1000):  
        is_reference = False
        for refseq_name in refseq_names:
            if refseq_name in cur_tag:
                is_reference = True
        sequences += [Seq(seq=cur_seq, name=cur_tag, is_reference=is_reference)]
    return sequences

# Build sequence logo
def get_sequence_logo(sequences):
    bcov_aln_counts = np.zeros((len(sequences[0].seq), 6))
    total_sequences = 0
    for sequence in sequences:
        if not sequence.is_shortlist:
            continue
        total_sequences += 1
        chars = list(sequence.seq)
        started = False # Avoid counting conservation until 5' end starts
        start_lastdash = 0 # Keeps track of number of continuous dashes at the end
        restart_lastdash = True
        for ii, cur_char in enumerate(chars):
            if cur_char != '-':
                started = True
                restart_lastdash = True
            elif restart_lastdash:
                start_lastdash = ii
                restart_lastdash = False
            cur_idx = 5
            if cur_char in pos_idx:
                cur_idx = pos_idx[cur_char]
            if started:
                bcov_aln_counts[ii, cur_idx] += 1
        if not restart_lastdash:
            for ii in range(start_lastdash, len(chars)):
                bcov_aln_counts[ii, 4] -= 1
    bcov_aln_counts = bcov_aln_counts/bcov_aln_counts.sum(axis=1, keepdims=True)
    return bcov_aln_counts

def get_aln_percs(aln_counts, reference_seq):
    aln_vals = []
    for ii, char in enumerate(reference_seq):
        if char != '-':
            aln_vals += [aln_counts[ii][pos_idx[char]]]
    return aln_vals

def get_intervals(aln_counts, sequences, cutoff=1, min_size=10):
    intervals = []
    start_int = -1
    end_int = -1
    in_int = False
    int_size = 0
    for ii in range(len(sequences[0].seq)):
        # Current working interval is over
        if max(aln_counts[ii, :]) < cutoff:
            if in_int:
                if int_size > min_size:
                    intervals += [(start_int, end_int, int_size)]
                int_size = 0
                in_int = False
        else:
            if not in_int:
                start_int = ii
                end_int = ii
                in_int = True
            if aln_counts[ii, 4] < cutoff:
                int_size += 1
            end_int = ii
    if in_int:
        if int_size > min_size:
            intervals += [(start_int, end_int, int_size)]
    return intervals

# Function to get positions and sequences in terms of the reference sequence, not the MSA
def get_seq(reference_seq):
    return reference_seq.replace('-', '')
        
def get_position(reference_seq, idx):
    reference_seq_l = list(reference_seq)
    cur_pos = 0
    for ii, char in enumerate(reference_seq_l):
        if ii == idx:
            break
        if char != '-':
            cur_pos += 1
    return cur_pos

def get_position_full(reference_seq, idx):
    reference_seq_l = list(reference_seq)
    cur_pos = 0
    for ii, char in enumerate(reference_seq_l):
        if cur_pos == idx:
            break
        if char != '-':
            cur_pos += 1
    return ii

def get_ref_seq(sequences):
    reference_seq = ''
    for sequence in sequences:
        if sequence.is_reference:
            reference_seq = sequence
    ref_seq = get_seq(reference_seq.seq) # Condensed (no gaps)
    return (reference_seq.seq, ref_seq)

def get_intervals_refseq(bcov_aln_counts, full_ref_seq, sequences, cutoff=1, MIN_SIZE=15):
    intervals = get_intervals(bcov_aln_counts, sequences, cutoff=cutoff, min_size=MIN_SIZE)
    ref_intervals = []
    for interval in intervals:
        start_pos = interval[0]
        end_pos = interval[1]

        ref_start = get_position(full_ref_seq, start_pos)
        ref_end = get_position(full_ref_seq, end_pos)
        if (ref_end - ref_start) >= MIN_SIZE:
            ref_intervals += [(ref_start + 1, ref_end + 1)]
        else:
            ref_intervals += [(-1, -1)]
    return ref_intervals

# Gets secondary structure
def get_secstruct_mea(int_start, int_end, ref_seq, secstruct_interval=20):
    sequence = ref_seq[(int_start-secstruct_interval):(int_end+secstruct_interval)]
    bp_matrix = bpps(sequence, package='contrafold_2')
    
    best_struct = ''
    best_mcc = 0
    for log_gamma in range(-2,2):
        mea_mdl = MEA(bp_matrix,gamma=10**log_gamma)
        [exp_sen, exp_ppv, exp_mcc, exp_fscore] = mea_mdl.score_expected()
        if exp_mcc > best_mcc:
            best_struct = mea_mdl.structure
            best_mcc = exp_mcc

    conserved_str = '.'*secstruct_interval + '*'*(int_end - int_start) + '.'*secstruct_interval
    return((best_struct, best_mcc, conserved_str, sequence))

# Generates Interval objects with secondary structure information for windows around each interval
def get_all_secstructs_mea(ref_intervals, ref_seq, secstruct_interval=20):
    all_intervals = []
    for ii in range(len(ref_intervals)):
        if ref_intervals[ii][0] == -1:
            continue
        int_start = ref_intervals[ii][0]
        int_end = ref_intervals[ii][1]
        (best_struct, best_mcc, conserved_str, full_seq) = get_secstruct_mea(int_start, int_end, ref_seq, secstruct_interval=secstruct_interval)
        new_int = Interval(seq=ref_seq[int_start:int_end], int_start=int_start, int_end=int_end, \
                secstruct=best_struct, mcc=best_mcc, conserved_region=conserved_str, full_seq=full_seq)
        all_intervals += [new_int]
    mcc_vals = np.array([x.mcc for x in all_intervals])
    # Sort in decreasing order of MCC
    sort_idxs = np.flip(np.argsort(mcc_vals))
    all_intervals = np.array(all_intervals)[sort_idxs]
    return all_intervals

# Get intervals based on reference sequence numbering from an alignment in filename that 
# have sequence conservation at least cutoff and size at least MIN_SIZE
def get_ref_intervals_from_file(filename, cutoff=1, MIN_SIZE=15):
    sequences = get_sequences(filename)
    bcov_aln_counts = get_sequence_logo(sequences)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)
    aln_vals = get_aln_percs(bcov_aln_counts, full_ref_seq)
    ref_intervals = get_intervals_refseq(bcov_aln_counts, full_ref_seq, sequences, cutoff=cutoff, MIN_SIZE=MIN_SIZE)
    return ref_intervals

def get_interval_overlap_single(test_interval, intervals): 
    has_overlap = False
    for interval in intervals:
        does_overlap = True
        if test_interval[0] < interval[0] and test_interval[1] < interval[0]:
            does_overlap = False
        if test_interval[0] > interval[1] and test_interval[1] > interval[1]:
            does_overlap = False
        if does_overlap:
            has_overlap = True
            break
    return has_overlap

# Get all intervals in intervals1 that overlap with some interval in intervals2
def get_interval_overlap(intervals1, intervals2):
    overlap_intervals = []
    
    for interval1 in intervals1:
        for interval2 in intervals2:
            does_overlap = True
            if (interval1[1] < interval2[0]) or (interval1[0] > interval2[1]):
                does_overlap = False
            if does_overlap:
                overlap_intervals += [interval1]
                break
    return overlap_intervals

# Get all intervals in intervals1 that overlap with some interval in intervals2
# by at least size nucleotides
def get_interval_overlap_size(intervals1, intervals2, min_size=1):
    overlap_intervals = []
    
    for interval1 in intervals1:
        max_overlap_size = 0
        for interval2 in intervals2:
            overlap_size = -1
            if (interval1[1] < interval2[0]) or (interval1[0] > interval2[1]):
                overlap_size = 0
            elif (interval1[1] >= interval2[0]) and (interval2[1] >= interval1[1]):
                overlap_size = min(interval1[1] - interval1[0] + 1, interval1[1] - interval2[0] + 1)
            elif (interval1[0] >= interval2[0]) and (interval1[0] <= interval2[1]):
                overlap_size = min(interval2[1] - interval1[0] + 1, interval1[1] - interval1[0] + 1)
            elif (interval1[0] <= interval2[0]) and (interval1[1] >= interval2[1]):
                overlap_size = interval2[1] - interval2[0]
            if overlap_size > max_overlap_size:
                max_overlap_size = overlap_size
            #if overlap_size > min_size:
            #    print(interval1)
            #    print(interval2)
        if max_overlap_size >= min_size:
            overlap_intervals += [interval1]
    return overlap_intervals

def get_num_overlaps_rnd_trials(intervals_list, intervals, refseq_len, num_trials=10000):
    lengths = [(interval[1] - interval[0]) for interval in intervals_list]
    num_overlaps = []
    for ii in range(num_trials):
        num_overlap = 0
        for length in lengths:
            start_idx = randint(0, refseq_len - length)
            if get_interval_overlap_single((start_idx, start_idx + length), intervals):
                num_overlap += 1
        num_overlaps += [num_overlap]
    return num_overlaps

def get_perc_overlaps_rnd_trials(intervals_list, intervals, refseq_len, num_trials=10000):
    lengths = [(interval[1] - interval[0]) for interval in intervals_list]
    perc_overlap = []
    for length in lengths:
        num_overlap = 0
        for ii in range(num_trials):
            start_idx = randint(0, refseq_len - length)
            if get_interval_overlap_single((start_idx, start_idx + length), intervals):
                num_overlap += 1
        perc_overlap += [num_overlap/num_trials]
    return perc_overlap

def get_num_overlaps_rnd_trials_size(regions_list, intervals, refseq_len, num_trials=10000, min_size=1):
    lengths = [(region[1] - region[0]) for region in regions_list]
    num_overlaps = []
    for ii in range(num_trials):
        num_overlap = 0
        for length in lengths:
            start_idx = randint(0, refseq_len - length)
            if get_interval_overlap_size([(start_idx, start_idx + length)], intervals, min_size=min_size):
                num_overlap += 1
        num_overlaps += [num_overlap]
    return num_overlaps

def convert_intervals_refseq(intervals, full_ref_seq):
    ref_intervals = []
    for interval in intervals:
        start_pos = interval[0]
        end_pos = interval[1]

        ref_start = get_position(full_ref_seq, start_pos)
        ref_end = get_position(full_ref_seq, end_pos)
        if len(interval) == 3: ## potentially included other interval data
            ref_intervals += [(ref_start + 1, ref_end + 1, interval[2], start_pos, end_pos)]
        else:
            ref_intervals += [(ref_start + 1, ref_end + 1, start_pos, end_pos)]
    return ref_intervals

# Write all RNAz windows and loci
def get_rnaz_windows(rnaz_windows_file, aln_file):
    f = open(rnaz_windows_file)
    rnaz_windows = f.readlines()
    f.close()

    rnaz_windows_dict = {}

    sequences = get_sequences(aln_file)
    (full_ref_seq, ref_seq) = get_ref_seq(sequences)

    ii = 0
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
            rnaz_windows_dict[(int_start, int_end)] = (seq, secstruct, z_score, p_val)
        ii += 1
    
    return rnaz_windows_dict