#!/usr/bin/env python
# -*- coding: utf-8 -*-

#IMPORTS AND FILE LOADING
import numpy as np
import random
import pandas as pd
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.Restriction import *
import json
from Bio.SeqUtils import gc_fraction
from tqdm import tqdm
import os
import argparse
import csv
import importlib.resources as pkg_resources

"""
This section loads the reference files needed for oligo assembly.
"""
ref_dir = pkg_resources.files('TCRAFT.references')
with open(ref_dir.joinpath('H_sapiens_codon_usage.json'), 'r') as infile:
    codon_freq_dict = json.load(infile)

trbv_trac_ohset = pd.read_csv(ref_dir.joinpath('TRBV_TRAC_OHset.csv'), index_col='N')
trbc_trav_ohset = pd.read_csv(ref_dir.joinpath('TRBC_TRAV_OHset.csv'), index_col='N')
trbv_trac_OH_frags = pd.read_csv(ref_dir.joinpath('TRBV_TRAC_OH_frags.csv'))
trbc_trav_OH_frags = pd.read_csv(ref_dir.joinpath('TRBC_TRAV_OH_frags.csv'))
trj_seqs = pd.read_csv(ref_dir.joinpath('TRJ_seqs.csv'))

with open(ref_dir.joinpath('Oligo_Components.json'), 'r') as infile:
    oligo_components = json.load(infile)

#Load primer and restriction site sequences
fwd_primer = oligo_components['fwd_primer']
rev_primer = oligo_components['rev_primer']
bbs1 = oligo_components['bbs1']
bsa1 = oligo_components['bsa1']
bsmb1 = oligo_components['bsmb1']


#DEFAULTS
#A fixed default seed makes runs reproducible: the same input CSV yields the same
#oligos, which makes results verifiable and failures replayable.
DEFAULT_SEED = 0


#Column names that are recognised automatically as holding clonotype names,
#in order of preference. The match is case-insensitive.
CLONOTYPE_CANDIDATE_COLUMNS = [
    'clonotype', 'clonotype_id', 'clonotype_name', 'clone_id', 'clone',
    'tcr_id', 'tcr_name', 'id', 'name', 'sample_id', 'label',
]

#Column created when the input CSV holds no clonotype name at all.
GENERATED_CLONOTYPE_COLUMN = 'clonotype'

#Required TCR annotation columns.
REQUIRED_COLUMNS = ['V_alpha', 'V_beta', 'CDR3_alpha', 'CDR3_beta', 'J_alpha', 'J_beta']


#HELPER FUNCTIONS
def detect_clonotype_column(columns):
    """
    Finds the column of an input CSV that holds clonotype names.

    params:
        columns: an iterable of column names, or a pandas DataFrame.

    return:
        The name of the detected column as it is spelled in the CSV, or None when
        no candidate column is present. Required TCR annotation columns are never
        returned.
    """
    if isinstance(columns, pd.DataFrame):
        columns = columns.columns

    columns = list(columns)
    lowered = {str(col).strip().lower(): col for col in columns}
    required = {col.lower() for col in REQUIRED_COLUMNS}

    for candidate in CLONOTYPE_CANDIDATE_COLUMNS:
        if candidate in lowered and candidate not in required:
            return lowered[candidate]

    return None


def list_input_columns(input_csv_path):
    """
    Returns the column names of an input CSV without loading the whole file.
    Used by the GUI to populate the clonotype column selector.
    """
    return list(pd.read_csv(input_csv_path, nrows=0).columns)


def generate_clonotype_names(count, width=4):
    """
    Builds default clonotype names (TCR_0001, TCR_0002, ...) for an input CSV
    that has no clonotype column, so that every oligo always carries a name.
    """
    width = max(width, len(str(count)))
    return [f'TCR_{i + 1:0{width}d}' for i in range(count)]


def summarize_input(input_csv_path):
    """
    Inspects an input CSV without generating anything, so the GUI can show the
    user what was understood from their file before a run is started.

    return:
        A dict with: n_rows, columns, clonotype_column (detected or None),
        missing_columns, n_duplicates, n_unique, error (None when readable).
    """
    summary = {'n_rows': 0, 'columns': [], 'clonotype_column': None,
               'missing_columns': list(REQUIRED_COLUMNS), 'n_duplicates': 0,
               'n_unique': 0, 'error': None}

    try:
        df = pd.read_csv(input_csv_path)
    except Exception as exc:
        summary['error'] = f'Could not read this CSV file: {exc}'
        return summary

    summary['n_rows'] = int(df.shape[0])
    summary['columns'] = [str(c) for c in df.columns]
    summary['clonotype_column'] = detect_clonotype_column(df)
    summary['missing_columns'] = [col for col in REQUIRED_COLUMNS if col not in df.columns]

    if summary['missing_columns'] or df.empty:
        return summary

    # Same normalisation as the run itself, so the preview matches the result.
    keys = pd.DataFrame({col: _normalize_tcr_field(df[col]) for col in REQUIRED_COLUMNS})
    summary['n_duplicates'] = int(keys.duplicated(keep='first').sum())
    summary['n_unique'] = summary['n_rows'] - summary['n_duplicates']

    return summary


def _uniquify_clonotype_names(series):
    """
    Makes clonotype names unique by suffixing repeats (name, name_2, name_3...).
    Two *different* TCRs sharing a name would otherwise produce two rows with the
    same label in the final table.

    return:
        A tuple (unique names as a Series, list of (original name, occurrences)).
    """
    names = series.astype(str)

    repeated = names.value_counts()
    repeated = repeated[repeated > 1]
    collisions = [(name, int(n)) for name, n in repeated.items()]

    if not collisions:
        return names, []

    # Suffixing can itself collide with an existing name, so repeat until stable.
    for _ in range(100):
        occurrence = names.groupby(names).cumcount()
        if not names.duplicated().any():
            break
        names = names.where(occurrence == 0, names + '_' + (occurrence + 1).astype(str))

    return names, collisions


def _normalize_tcr_field(series):
    """
    Normalises a V/J/CDR3 column for duplicate detection: strips surrounding and
    inner whitespace and upper-cases the value. Used only to build the duplicate
    key; the original values are preserved in all output files.
    """
    return series.astype(str).str.replace(r'\s+', '', regex=True).str.upper()


def rev_comp(seq):
    """
    Returns the reverse complement of a DNA sequence 'seq'.
    """
    mapper = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    return ''.join([mapper[nt] for nt in seq[::-1]])

def dfs_seq_builder(ref_aa, primer='', aa_pos=0):
    """
    This function uses a depth-first search to build all possible DNA sequences that translate to the given amino acid sequence.

    params:
        ref_aa: str, the reference amino acid sequence
        primer: str, the current DNA sequence being built
        aa_pos: int, the current position in the amino acid sequence

    returns:
        list of all possible DNA sequences that translate to the given amino acid sequence
    """
    if aa_pos >= len(ref_aa):
        return [primer]
    output = []
    for codon_option, _ in codon_freq_dict[ref_aa[aa_pos]]:
        output.extend(dfs_seq_builder(ref_aa, primer+codon_option, aa_pos+1))
    return output 


#FUNCTIONS FOR CDR3 CONSTRUCTION
def check_re_sites(seq, sites):
    """
    Checks if any of the given restriction sites are present in the sequence.
    params
        seq: The sequence to check.
        sites: A list of restriction sites to check for.

    return
        A tuple containing the position of the first found restriction site and the site itself.
        If no sites are found, returns (-1, None).
    """
    for site in sites:
        if site in seq or rev_comp(site) in seq:
            try:
                pos = seq.index(site)
            except ValueError:
                pos = seq.index(rev_comp(site))
            return pos, site
    return -1, None


def check_homopolymer_sites(seq):
    """
    Checks if homopolymers >= 5bp in length are present in the sequence.
    params
        seq: The sequence to check.

    return
        A tuple containing the position of the first found homopolymer site and the site itself.
        If no sites are found, returns (-1, None).
    """
    for site in ['AAAAA', 'CCCCC', 'GGGGG', 'TTTTT']:
        if site in seq:
            try:
                pos = seq.index(site)
            except ValueError:
                pos = seq.index(rev_comp(site))
            return pos, site
    return -1, None


def remove_restriction_sites(seq, aa, sites):
    """
    Iteratively removes restriction sites from the sequence by replacing them with a random codon that does not contain the restriction site.
    params:
        seq: The sequence to check.
        aa: The amino acid sequence corresponding to the DNA sequence.
        sites: A list of restriction sites to check for and remove.

    return:
        The modified sequence with restriction sites removed.
    """
    iters = 0
    while check_re_sites(seq, sites)[0] != -1:
        pos, site = check_re_sites(seq, sites)
        start = max(0, pos // 3 - 2)
        end = min((pos + len(site)) // 3 + 2, len(seq) // 3 - 1)       
         
        replace = None
        for option in dfs_seq_builder(aa[start:end+1]):
            if check_re_sites(option, sites)[0] == -1:
                replace = option

        if iters >= 15 or replace is None:
            raise ValueError(f'Note: Failed to remove restriction site. Seq = {seq}.\nThis usually happens due to a bad random seed choice and can often be resolved by rerunning the script.')
        
        seq = seq[:start*3] + replace + seq[(end+1)*3:]  
        iters += 1
    
    assert Seq(seq).translate() == aa
    return seq


def remove_homopolymer_sites(seq, aa, sites):
    """
    Iteratively removes homopolymer sites from the sequence by replacing them with a random codon that does not contain the homopolymer site.
    Simultaneously checks for restriction sites to ensure they are not reintroduced after removing homopolymers.

    params:
        seq: The sequence to check.
        aa: The amino acid sequence corresponding to the DNA sequence.
        sites: A list of restriction sites to check for and remove. This is to ensure that the homopolymer replacement does not introduce new restriction sites.

    return:
        The modified sequence with restriction sites removed.
    """
    iters = 0
    while check_homopolymer_sites(seq)[0] != -1:
        pos, site = check_homopolymer_sites(seq)
        start = max(0, pos // 3 - 2)
        end = min((pos + len(site)) // 3 + 2, len(seq) // 3 - 1)
        
        replace = None
        for option in dfs_seq_builder(aa[start:end+1]):
            if check_re_sites(option, sites)[0] == -1 and check_homopolymer_sites(option)[0] == -1:
                replace = option 

        if iters >= 15 or replace is None:
            raise ValueError(f'Note: Failed to remove homopolymer. Seq = {seq}\nnThis usually happens due to a bad random seed choice and can often be resolved by rerunning the script.')
        
        seq = seq[:start*3] + replace + seq[(end+1)*3:]
        iters += 1
        
    assert Seq(seq).translate() == aa
    return seq


def generate_cleaned_DNA(aa_seq, freq_cutoff=0.15):
    """
    Generates a cleaned DNA sequence from the given amino acid sequence.
    The DNA sequence is generated by sampling from the codon usage table, ensuring that the GC content is between 25% and 65%.
    The generated sequence is then cleaned by removing restriction sites and homopolymer sites.
    params:
        aa_seq: str, the amino acid sequence to generate DNA for.
        freq_cutoff: float, the minimum frequency of a codon to be considered for sampling.

    return:
        The cleaned DNA sequence as a string
    """
    
    def generate_DNA(): #helper function to generate a DNA sequence for the given amino acid sequence.
        seq = ""
        while seq == "" or gc_fraction(seq) < 0.25 or gc_fraction(seq) > 0.65:
            seq = ""
            for aa in aa_seq:
                codons = []
                freqs = []
                for codon, freq in codon_freq_dict[aa]:
                    if freq > freq_cutoff:
                        codons.append(codon)
                        freqs.append(freq)
                freqs = (np.array(freqs) / np.sum(freqs)).tolist()
                seq += random.choices(population=codons, weights=freqs, k=1)[0]
        assert Seq(seq).translate() == aa_seq
        return seq

    #clean the DNA
    tries = 0
    while True:
        try:
            seq = generate_DNA()
            sites = [BsaI.site, BsmBI.site, BbsI.site, SapI.site, XbaI.site, EcoRI.site, XhoI.site, EcoRV.site, NotI.site]
            seq_minus_re = remove_restriction_sites(seq, aa_seq, sites)
            final_seq = remove_homopolymer_sites(seq_minus_re, aa_seq, sites)
            assert check_re_sites(final_seq, sites)[0] == -1
            assert check_homopolymer_sites(final_seq)[0] == -1
            break
        except Exception as e:
            tries += 1
            if tries > 25:
                print(e)
                raise ValueError('Note: Failed to generate cleaned DNA after 25 sampling attempts.\nThis usually happens due to a bad random seed choice and can often be resolved by rerunning the script.')
    
    return final_seq


def concatenate_oligo(oligo):
    """
    Concatenates the oligo fragments into a single sequence.
    params:
        oligo: dict, the oligo fragments to concatenate.
    return:
        The concatenated oligo sequence as a string.
    """
    final_oligo_seq = fwd_primer + \
                      oligo['left_re'] + \
                      oligo['trbv_right_frag'] + \
                      oligo['cdr3bj_seq'] + \
                      oligo['trbc_left_frag'] + \
                      oligo['middle_re'] + \
                      oligo['trav_right_frag'] + \
                      oligo['cdr3aj_seq'] + \
                      oligo['trac_left_frag'] + \
                      oligo['right_re'] + \
                      rev_primer
    return final_oligo_seq


def test_oligo(oligo_seq, trav_set):
    """
    Tests the fully assembled oligo sequence for restriction sites.
    params:
        oligo_seq: str, the oligo sequence to test.
        trav_set: str, the TRAV set to check for the right restriction sites.
    return:
        True if the oligo sequence passes the tests, False otherwise.
    """
    # check restriction sites
    rb = BsmBI + BsaI + BbsI + SapI + XbaI + EcoRI + XhoI + NotI + EcoRV
    restriction_analysis = Analysis(rb, Seq(oligo_seq)).full()
    
    try:
        for key in restriction_analysis.keys():
            if key == BbsI: #need to patch a structural problem with trav12-3
                assert len(restriction_analysis[key]) == 2
            elif key == BsaI:
                if trav_set == 'A':
                    assert len(restriction_analysis[key]) == 1
                else: 
                    assert len(restriction_analysis[key]) == 2
            elif key == BsmBI:
                if trav_set == 'A':
                    assert len(restriction_analysis[key]) == 2
                else: 
                    assert len(restriction_analysis[key]) == 1
            else: #all others minus XhoI, XbaI
                assert len(restriction_analysis[key]) == 0
    except AssertionError:
        return False
    return True


def create_CDR3_oligo(trav_name, trbv_name, cdr3a, cdr3b, traj_name, trbj_name):
    """
    Creates a CDR3 oligo from the given TCR information.
    params:
        trav_name: TRAV allele ID
        trbv_name: TRBV allele ID
        cdr3a: CDR3 alpha amino acid sequence, starting from C and ending with F
        cdr3b: CDR3 beta amino acid sequence, starting from C and ending with F
        traj_name: TRAJ allele ID
        trbj_name: TRBJ allele ID
    return:
        A tuple containing fully assembled oligo sequence as a string and the TRBV set for Step 1 pooling
    """
    oligo_frags = {}

    if trav_name not in trbc_trav_OH_frags['trav'].values or trbv_name not in trbv_trac_OH_frags['trbv'].values:
        raise ValueError(f'Note: Invalid TRAV {trav_name} and/or TRBV {trbv_name}, skipping this TCR.')
    
    oligo_frags['trbv_right_frag'] = trbv_trac_OH_frags['trbv_right_frag'].loc[trbv_trac_OH_frags['trbv'] == trbv_name].values[0]
    oligo_frags['trbc_left_frag'] = trbc_trav_OH_frags['trbc_left_frag'].loc[trbc_trav_OH_frags['trav'] == trav_name].values[0]
    oligo_frags['trav_right_frag'] = trbc_trav_OH_frags['trav_right_frag'].loc[trbc_trav_OH_frags['trav'] == trav_name].values[0]
    oligo_frags['trac_left_frag'] = trbv_trac_OH_frags['trac_left_frag'].loc[trbv_trac_OH_frags['trbv'] == trbv_name].values[0]

    #Find TRAJ and TRBJ
    try:
        traj = trj_seqs.loc[trj_seqs.TRJ == traj_name, 'seq'].values[0]
        trbj = trj_seqs.loc[trj_seqs.TRJ == trbj_name, 'seq'].values[0]

    except IndexError:
        raise ValueError(f'Note: Invalid TRAJ {traj_name} and/or TRBJ {trbj_name}, skipping this TCR.')
    

    #prepare full CDR3A+J and CDR3B+J seq, cut off the first C (already in the TRAV/TRBV)
    #residue between TRAJ and TRAC is specified in the TRAJ
    oligo_frags['cdr3aj_seq'] = generate_cleaned_DNA(cdr3a[1:] + traj).lower()
    oligo_frags['cdr3bj_seq'] = generate_cleaned_DNA(cdr3b[1:] + trbj).lower()

    #patch frequent edge cases that occur in the J-C junction to prevent excessive resampling
    if oligo_frags['cdr3bj_seq'][-3:] == 'ctc': #patch frequent CTCGAG edge case in TRBJ-TRBC junction
        oligo_frags['cdr3bj_seq'] = oligo_frags['cdr3bj_seq'][:-3] + 'ctg' #another leucine codon 
    if oligo_frags['cdr3aj_seq'][-3:] == 'gat': #patch frequent GATATC edge case in TRAJ-TRAC junction
        oligo_frags['cdr3aj_seq'] = oligo_frags['cdr3aj_seq'][:-3] + 'gac' #another aspartate codon

    #generate restriction site sequences
    trbv_set = trbv_trac_ohset.loc[trbv_trac_ohset.TRBV_name == trbv_name, 'set'].values[0]
    trav_set = trbc_trav_ohset.loc[trbc_trav_ohset.TRAV_name == trav_name, 'set'].values[0]

    # all CDRs have BbsI for step 1
    oligo_frags['left_re'], oligo_frags['right_re'] = bbs1

    # assign middle restriction site based on TRAV set
    if trav_set == 'A':
        oligo_frags['middle_re'] = bsmb1
    else:
        oligo_frags['middle_re'] = bsa1

    #resample CDR3s if restriction sites are present in the junctions between sequence fragments
    resample_count = 0
    while not test_oligo(concatenate_oligo(oligo_frags), trav_set):
        oligo_frags['cdr3aj_seq'] = generate_cleaned_DNA(cdr3a[1:] + traj).lower()
        oligo_frags['cdr3bj_seq'] = generate_cleaned_DNA(cdr3b[1:] + trbj).lower()
        resample_count += 1
        if resample_count > 25:
            raise ValueError('Note: Failed to fix restriction sites in oligo junctions after 25 resampling attempts. Skipping this TCR for now.\n ' \
            'This usually happens due to a bad random seed choice and can often be resolved by rerunning the script. ')
        
    if resample_count > 0:
        print(f'Resampled {resample_count} times to fix restriction site issue.')

    return concatenate_oligo(oligo_frags), trbv_set #to split Step1 pools


# PADDING FUNCTIONS
TYPE_IIS_PADDING_SITES={"SapI":["GCTCTTC","GAAGAGC"],"BsaI":["GGTCTC","GAGACC"],"BbsI":["GAAGAC","GTCTTC"],"BsmBI":["CGTCTC","GAGACG"]}

def _find_sites(seq,site):
    seq,site=seq.upper(),site.upper(); positions=[]; start=0
    while True:
        pos=seq.find(site,start)
        if pos<0: return positions
        positions.append(pos); start=pos+1

def _scan_padding_sites(seq):
    result={}
    for enzyme,sites in TYPE_IIS_PADDING_SITES.items():
        for site in sites:
            positions=_find_sites(seq,site)
            if positions: result.setdefault(enzyme,{})[site]=positions
    return result

def _find_homopolymers(seq,min_length=5):
    result=[]; seq=seq.upper()
    for base in "ACGT":
        start=0; site=base*min_length
        while True:
            pos=seq.find(site,start)
            if pos<0: break
            result.append((base,pos)); start=pos+1
    return result

def _new_padding_sites(original,padded,left_len):
    new=[]; original_len=len(original)
    for enzyme,sites in _scan_padding_sites(padded).items():
        for site,positions in sites.items():
            original_positions=set(_find_sites(original,site)); site_len=len(site)
            for pos in positions:
                original_pos=pos-left_len
                if original_pos>=0 and original_pos+site_len<=original_len:
                    if original_pos not in original_positions: new.append((enzyme,site,pos,"original_region"))
                else:
                    location="left_padding_or_junction" if pos<left_len else "right_padding_or_junction" if pos>=left_len+original_len else "junction"
                    new.append((enzyme,site,pos,location))
    return new

def _new_padding_homopolymers(original,padded,left_len):
    new=[]; original_len=len(original); original_positions={(b,p) for b,p in _find_homopolymers(original)}
    for base,pos in _find_homopolymers(padded):
        original_pos=pos-left_len
        if original_pos>=0 and original_pos+5<=original_len:
            if (base,original_pos) not in original_positions: new.append((base,pos,"original_region"))
        else:
            location="left_padding_or_junction" if pos<left_len else "right_padding_or_junction" if pos>=left_len+original_len else "junction"
            new.append((base,pos,location))
    return new

def _find_safe_padding(original, target, random_seq, preferred_offset=0):
    """
    Finds a 3' padding sequence that brings 'original' up to 'target' nt without
    creating a new Type-IIS restriction site or a homopolymer (>=5 nt), neither
    inside the padding itself nor at the junction with the end of the oligo.

    The padding is always taken from the single shared random reference sequence,
    so every oligo of a pool receives the SAME padding whenever that is safe
    (preferred_offset, normally 0). Only when that shared padding would create a
    site or a homopolymer for this particular oligo is the window slid along the
    reference sequence until a safe one is found.

    params:
        original: str, the unpadded oligo sequence.
        target: int, the desired final length in nt.
        random_seq: str, the shared random reference sequence (ACGT only).
        preferred_offset: int, the window offset to try first.

    return:
        A dict describing the padding, or None when no safe window exists.
    """
    needed = target - len(original)

    if needed < 0:
        raise ValueError(f'Oligo is {-needed} nt longer than the padding target of {target} nt.')

    if needed == 0:
        return {'padding': '', 'padded': original, 'offset': 0,
                'needed': 0, 'status': 'not_needed'}

    if needed > len(random_seq):
        raise ValueError(f"Need {needed} nt of padding but random sequence contains only {len(random_seq)} nt.")

    max_offset = len(random_seq) - needed
    preferred_offset = min(max(preferred_offset, 0), max_offset)

    offsets = [preferred_offset] + [o for o in range(max_offset + 1) if o != preferred_offset]

    for offset in offsets:
        padding = random_seq[offset:offset + needed]
        padded = original + padding

        if _new_padding_sites(original, padded, 0):
            continue
        if _new_padding_homopolymers(original, padded, 0):
            continue

        return {'padding': padding, 'padded': padded, 'offset': offset, 'needed': needed,
                'status': 'padded' if offset == preferred_offset else 'padded_shifted_window'}

    return None


def _pad_row_label(row, clonotype_column=None):
    """
    Builds a human-readable label for a row of a pool table, used in the padding report.
    """
    if clonotype_column is not None and clonotype_column in row and pd.notna(row[clonotype_column]):
        return str(row[clonotype_column])

    detected = detect_clonotype_column(list(row.index))
    if detected is not None and pd.notna(row[detected]):
        return str(row[detected])

    return f"{row.get('TRAV', '?')}|{row.get('CDR3A', '?')}|{row.get('TRBV', '?')}|{row.get('CDR3B', '?')}"


def _pad_pool_table(df, random_seq, table_name='', clonotype_column=None):
    """
    Pads every oligo of a pool table up to the longest oligo of that same table.
    Oligos for which no safe padding exists are left unpadded and reported; the
    run is never aborted because of a single oligo.

    return:
        A tuple (padded dataframe, list of report rows).
    """
    df = df.copy()

    if df.empty:
        df['Padding_added'] = pd.Series(dtype='int64')
        df['Padding_status'] = pd.Series(dtype='object')
        return df, []

    target = int(df['Sequence'].str.len().max())

    added = []
    statuses = []
    report = []

    for _, row in df.iterrows():
        original = str(row['Sequence'])
        result = _find_safe_padding(original, target, random_seq)

        if result is None:
            added.append(0)
            statuses.append('failed_no_safe_padding')
            report.append({'Table': table_name, 'Oligo': _pad_row_label(row, clonotype_column),
                           'Original_length': len(original), 'Target_length': target,
                           'Padding_needed': target - len(original), 'Padding_added': 0,
                           'Status': 'failed_no_safe_padding'})
            continue

        added.append(result['needed'])
        statuses.append(result['status'])

        if result['status'] != 'not_needed':
            df.loc[row.name, 'Sequence'] = result['padded']
            df.loc[row.name, 'Length'] = len(result['padded'])
            df.loc[row.name, 'GC'] = gc_fraction(result['padded'])

        if result['status'] == 'padded_shifted_window':
            report.append({'Table': table_name, 'Oligo': _pad_row_label(row, clonotype_column),
                           'Original_length': len(original), 'Target_length': target,
                           'Padding_needed': result['needed'], 'Padding_added': result['needed'],
                           'Status': f"padded_shifted_window (offset {result['offset']})"})

    df['Padding_added'] = added
    df['Padding_status'] = statuses

    return df, report


def _save_padded_tables(output_dir, random_sequence_file, clonotype_column=None):
    """
    Pads each pool table independently and writes the *_padded.csv files plus a
    Padding_report.csv listing every oligo that needed a shifted window or that
    could not be padded safely.
    """
    random_sequence_file = str(random_sequence_file)

    if not os.path.exists(random_sequence_file):
        raise FileNotFoundError(f'Random sequence file not found: {random_sequence_file}')

    with open(random_sequence_file, 'r') as f:
        random_seq = "".join(f.read().split()).upper()

    if set(random_seq) - set("ACGT"):
        raise ValueError(f'Invalid characters in random sequence: {set(random_seq)-set("ACGT")}')

    files = ['Under350_PoolA_CDR3_oligos.csv', 'Over350_PoolA_CDR3_oligos.csv',
             'Under350_PoolB_CDR3_oligos.csv', 'Over350_PoolB_CDR3_oligos.csv']

    full_report = []
    failed_total = 0

    for filename in files:
        path = os.path.join(output_dir, filename)

        if not os.path.exists(path):
            continue

        df = pd.read_csv(path)

        if df.empty:
            continue

        padded, report = _pad_pool_table(df, random_seq, table_name=filename,
                                         clonotype_column=clonotype_column)
        padded.to_csv(os.path.join(output_dir, filename.replace('.csv', '_padded.csv')), index=False)

        full_report.extend(report)
        failed = int((padded['Padding_status'] == 'failed_no_safe_padding').sum())
        shifted = int((padded['Padding_status'] == 'padded_shifted_window').sum())
        failed_total += failed

        print(f'Padded {filename}: target={int(padded["Sequence"].str.len().max())} bp, '
              f'{shifted} shifted window(s), {failed} not padded')

    report_df = pd.DataFrame(full_report, columns=['Table', 'Oligo', 'Original_length', 'Target_length',
                                                   'Padding_needed', 'Padding_added', 'Status'])
    report_df.to_csv(os.path.join(output_dir, 'Padding_report.csv'), index=False)

    if failed_total > 0:
        print(f'Warning: {failed_total} oligo(s) could not be padded safely and were left unpadded. '
              f'See Padding_report.csv.')


def generate_cdr3_oligos(input_csv_path, output_dir=None, clonotype_id_column=None,
                         random_sequence_file=None, seed=DEFAULT_SEED,
                         progress_callback=None):
    """
    Main function to generate CDR3 oligos from TCR information.
    params:
        input_csv_path: str, path to the input CSV file containing TCR information.
        output_dir: str, path to the output directory to save generated oligos.
                    Default is './TCRAFT-generate_<date>'
        clonotype_id_column: str, input CSV column holding clonotype names. When
                    None, a known column name is auto-detected; if the CSV holds
                    none, names TCR_0001... are generated. Every oligo therefore
                    always carries a clonotype name in every output table.
        random_sequence_file: path to the shared random sequence used for 3' padding.
                    Default is the bundled TCRAFT/references/random_1000nt.txt.
        seed: int or None, random seed. A fixed seed makes a run reproducible:
                    the same input CSV yields the same oligos. Pass None for a
                    non-reproducible run.
        progress_callback: optional callable(done, total) invoked after each TCR.
                    Used by the GUI to show real progress; when given, the tqdm
                    progress bar is suppressed.

    return:
        None, outputs files to the output directory.
    """

    print('*******************')
    print('* TCRAFT-generate *')
    print('*******************')
    print()

    #PREPARE INPUTS
    if not os.path.exists(input_csv_path):
        raise FileNotFoundError(f'Input CSV file {input_csv_path} not found.')

    if not input_csv_path.endswith('.csv'):
        raise ValueError('Input file must be a CSV file.')

    # The padding reference ships with the package: resolve it here so that the
    # GUI and direct imports work exactly like the CLI.
    if random_sequence_file is None:
        random_sequence_file = ref_dir.joinpath('random_1000nt.txt')

    # A fixed seed makes codon sampling reproducible across runs.
    if seed is not None:
        random.seed(seed)
        np.random.seed(int(seed) % (2 ** 32))
        print(f'Random seed: {seed}')

    if output_dir is None:
        date_str = pd.Timestamp.now().strftime('%Y%m%d')
        output_dir = os.path.join(os.getcwd(), f"TCRAFT-generate_{date_str}")

    # Create the output directory before anything is written into it.
    os.makedirs(output_dir, exist_ok=True)

    print('Loading TCRs and deduplicating...')

    tcr_list = pd.read_csv(input_csv_path)
    required_columns = list(REQUIRED_COLUMNS)

    # Check the input columns BEFORE using them to deduplicate, so a malformed
    # CSV produces a clear message instead of a KeyError.
    missing_columns = [col for col in required_columns if col not in tcr_list.columns]
    if missing_columns:
        raise ValueError(f'Input CSV is missing the following required column(s): {", ".join(missing_columns)}. '
                         f'Required columns are: {", ".join(required_columns)}')

    # RESOLVE THE CLONOTYPE NAME COLUMN
    # Explicit choice wins; otherwise a known column name is detected; otherwise
    # names are generated, so that every oligo always carries a clonotype name.
    if clonotype_id_column is not None:
        if clonotype_id_column not in tcr_list.columns:
            raise ValueError(f'Clonotype ID column "{clonotype_id_column}" was not found in the input CSV. '
                             f'Available columns are: {", ".join(str(c) for c in tcr_list.columns)}')
        print(f'Clonotype names: using column "{clonotype_id_column}".')

    else:
        detected = detect_clonotype_column(tcr_list)

        if detected is not None:
            clonotype_id_column = detected
            print(f'Clonotype names: auto-detected column "{clonotype_id_column}".')

        else:
            clonotype_id_column = GENERATED_CLONOTYPE_COLUMN
            tcr_list.insert(0, clonotype_id_column, generate_clonotype_names(tcr_list.shape[0]))
            print(f'Clonotype names: no clonotype column found in the input CSV, '
                  f'generated names TCR_0001...TCR_{tcr_list.shape[0]:04d} in column "{clonotype_id_column}".')

    # Rows with no name at all would otherwise end up as 'nan' in the output.
    unnamed = tcr_list[clonotype_id_column].isna()
    if unnamed.any():
        fallback = generate_clonotype_names(tcr_list.shape[0])
        tcr_list.loc[unnamed, clonotype_id_column] = [fallback[i] for i in np.flatnonzero(unnamed)]
        print(f'Warning: {int(unnamed.sum())} row(s) had an empty clonotype name; '
              f'generated TCR_#### names for them.')

    # Normalised copies of the VDJ columns (trimmed, upper-cased, inner whitespace
    # removed) are used as the duplicate key, so that 'TRAV6 ' and 'trav6' are
    # recognised as the same clonotype. The original values are left untouched in
    # every output file.
    key_columns = []
    for col in required_columns:
        key_col = f'__key_{col}'
        tcr_list[key_col] = _normalize_tcr_field(tcr_list[col])
        key_columns.append(key_col)

    # VDJ duplicate groups
    dup = tcr_list.duplicated(key_columns, keep=False)
    grp = tcr_list.groupby(key_columns, dropna=False).ngroup() + 1
    first = ~tcr_list.duplicated(key_columns, keep="first")

    # Status + kept clonotype
    tcr_list["dedup_status"] = "unique"
    tcr_list.loc[dup & ~first, "dedup_status"] = "removed"
    tcr_list.loc[dup & first, "dedup_status"] = "kept"
    tcr_list["VDJ_group"] = grp

    # Every removed duplicate records the name of the clonotype that was kept in
    # its place, so the audit file stays traceable.
    tcr_list["kept_clonotype"] = (tcr_list[clonotype_id_column].where(first).groupby(grp).transform("first"))

    # Outputs
    helper_columns = key_columns + ["dedup_status", "VDJ_group", "kept_clonotype"]
    tcr_list_kept = tcr_list[first].drop(columns=helper_columns, errors="ignore")
    tcr_list_audit = (tcr_list[dup].assign(action=lambda x: x["dedup_status"].str.upper()).sort_values(["VDJ_group", "dedup_status"]).drop(columns=key_columns + ["dedup_status"], errors="ignore"))

    removed_num = int((dup & ~first).sum())
    print(f'Loaded {tcr_list.shape[0]} TCRs, kept {tcr_list_kept.shape[0]} unique, '
          f'removed {removed_num} duplicate(s).')

    tcr_list = tcr_list_kept.copy()

    # Two DIFFERENT TCRs may still share the same clonotype name. Make the names
    # unique so that one name means one oligo in the final tables.
    unique_names, name_collisions = _uniquify_clonotype_names(tcr_list[clonotype_id_column])
    tcr_list[clonotype_id_column] = unique_names

    if name_collisions:
        print(f'Warning: {len(name_collisions)} clonotype name(s) were used by several distinct TCRs '
              f'and have been suffixed to stay unique:')
        for name, occurrences in name_collisions[:10]:
            print(f'  "{name}" appeared {occurrences} times -> "{name}", "{name}_2", ...')
        if len(name_collisions) > 10:
            print(f'  ... and {len(name_collisions) - 10} more (see All_TCR_deduplicated.csv).')

    # Write the deduplicated list after the names have been made unique, so the
    # saved table matches the names used in the oligo tables.
    tcr_list.to_csv(os.path.join(output_dir, f'All_TCR_deduplicated.csv'), index=False)
    tcr_list_audit.to_csv(os.path.join(output_dir, f'All_TCR_removed.csv'), index=False)

    input_columns = [clonotype_id_column] + required_columns
    tcr_list = tcr_list[input_columns].dropna()

    print('Generating CDR3 oligos...')
    
    #MAIN LOOP - GENERATE CDR3 OLIGOS
    invalid_cdr3_num = 0
    cdr3_oligo_list = {'TRAV':[], 'TRBV':[], 'CDR3A':[], 'CDR3B':[], 'TRAJ':[], 'TRBJ':[], 'Length':[], 'GC':[], 'Pool':[], 'Sequence':[]}
    invalid_cdr3_oligo_list = {'TRAV':[], 'TRBV':[], 'CDR3A':[], 'CDR3B':[], 'TRAJ':[], 'TRBJ':[], 'Issue':[]}
    if clonotype_id_column is not None:
        cdr3_oligo_list={clonotype_id_column:[]}|cdr3_oligo_list
        invalid_cdr3_oligo_list={clonotype_id_column:[]}|invalid_cdr3_oligo_list

    total_tcrs = tcr_list.shape[0]

    for i in tqdm(range(total_tcrs), disable=progress_callback is not None):
        tcr = tcr_list.iloc[i, :]
        try:
            cdr3_oligo, pool = create_CDR3_oligo(tcr['V_alpha'], tcr['V_beta'], tcr['CDR3_alpha'], tcr['CDR3_beta'], tcr['J_alpha'], tcr['J_beta'])
            if clonotype_id_column is not None:
                cdr3_oligo_list[clonotype_id_column].append(tcr[clonotype_id_column])
        
            cdr3_oligo_list['TRAV'].append(tcr['V_alpha'])
            cdr3_oligo_list['TRBV'].append(tcr['V_beta'])
            cdr3_oligo_list['CDR3A'].append(tcr['CDR3_alpha'])
            cdr3_oligo_list['CDR3B'].append(tcr['CDR3_beta'])
            cdr3_oligo_list['TRAJ'].append(tcr['J_alpha'])
            cdr3_oligo_list['TRBJ'].append(tcr['J_beta'])
            cdr3_oligo_list['Length'].append(len(cdr3_oligo))
            cdr3_oligo_list['GC'].append(gc_fraction(cdr3_oligo))
            cdr3_oligo_list['Pool'].append(pool)
            cdr3_oligo_list['Sequence'].append(cdr3_oligo)

        except ValueError as e:
            print(e)
            invalid_cdr3_num += 1
            
            if clonotype_id_column is not None:
                invalid_cdr3_oligo_list[clonotype_id_column].append(tcr[clonotype_id_column])

            invalid_cdr3_oligo_list['TRAV'].append(tcr['V_alpha'])
            invalid_cdr3_oligo_list['TRBV'].append(tcr['V_beta'])
            invalid_cdr3_oligo_list['CDR3A'].append(tcr['CDR3_alpha'])
            invalid_cdr3_oligo_list['CDR3B'].append(tcr['CDR3_beta'])
            invalid_cdr3_oligo_list['TRAJ'].append(tcr['J_alpha'])
            invalid_cdr3_oligo_list['TRBJ'].append(tcr['J_beta'])
            invalid_cdr3_oligo_list['Issue'].append(str(e))

        if progress_callback is not None:
            progress_callback(i + 1, total_tcrs)

    print(f'Generated {len(cdr3_oligo_list["Sequence"])} valid CDR3 oligos.')
    if invalid_cdr3_num > 0:
        print(f'Failed to generate {invalid_cdr3_num} CDR3 oligos, see log for details.')

    #POST-PROCESSING
    cdr3_oligo_list = pd.DataFrame(cdr3_oligo_list)
    invalid_cdr3_oligo_list = pd.DataFrame(invalid_cdr3_oligo_list)
    if cdr3_oligo_list.empty:
        invalid_cdr3_oligo_list.to_csv(os.path.join(output_dir, 'Invalid_CDR3_table.csv'), index=False)
        raise ValueError('No valid CDR3 oligo could be generated from this input. '
                         f'See Invalid_CDR3_table.csv in {output_dir} for the reason of each failure.')

    cdr3_oligo_list.sort_values('Pool', inplace=True)
    over_350_num = int(np.sum(cdr3_oligo_list['Length'] > 350))
    print(f'{over_350_num}/{cdr3_oligo_list.shape[0]} oligos ({round(100*over_350_num/cdr3_oligo_list.shape[0], 2)}%) are over 350 bp long.')

    poolA_under350 = cdr3_oligo_list.loc[np.logical_and(cdr3_oligo_list['Pool'] == 'A', cdr3_oligo_list['Length'] <= 350)].copy()
    poolA_over350 = cdr3_oligo_list.loc[np.logical_and(cdr3_oligo_list['Pool'] == 'A', cdr3_oligo_list['Length'] > 350)].copy()
    poolB_under350 = cdr3_oligo_list.loc[np.logical_and(cdr3_oligo_list['Pool'] == 'B', cdr3_oligo_list['Length'] <= 350)].copy()
    poolB_over350 = cdr3_oligo_list.loc[np.logical_and(cdr3_oligo_list['Pool'] == 'B', cdr3_oligo_list['Length'] > 350)].copy()
    all_under350 = cdr3_oligo_list.loc[cdr3_oligo_list['Length'] <= 350].copy()
    all_over350 = cdr3_oligo_list.loc[cdr3_oligo_list['Length'] > 350].copy()


    #SAVE OUTPUTS
    print(f'Saving output files to {output_dir}')
    os.makedirs(output_dir, exist_ok=True)


    cdr3_oligo_list.to_csv(os.path.join(output_dir, f'All_{cdr3_oligo_list.shape[0]}_CDR3_oligos.csv'), index=False)
    invalid_cdr3_oligo_list.to_csv(os.path.join(output_dir, 'Invalid_CDR3_table.csv'), index=False)
    
    all_under350.to_csv(os.path.join(output_dir, f'Under350_All_CDR3_oligos.csv'), index=False)
    all_over350.to_csv(os.path.join(output_dir, f'Over350_All_CDR3_oligos.csv'), index=False)

    poolA_under350.to_csv(os.path.join(output_dir, 'Under350_PoolA_CDR3_oligos.csv'), index=False)
    poolB_under350.to_csv(os.path.join(output_dir, 'Under350_PoolB_CDR3_oligos.csv'), index=False)
    poolA_over350.to_csv(os.path.join(output_dir, 'Over350_PoolA_CDR3_oligos.csv'), index=False)
    poolB_over350.to_csv(os.path.join(output_dir, 'Over350_PoolB_CDR3_oligos.csv'), index=False)

    # Pad each Pool A/B and Under/Over350 table independently.
    _save_padded_tables(output_dir, random_sequence_file, clonotype_column=clonotype_id_column)

    #save oligo pool metadata
    with open(os.path.join(output_dir, 'Oligo_Pool_Metadata.csv'), 'w') as out:
        out.write('Entry Name,Value\n')

        out.write(f'Random Seed,{seed if seed is not None else "none (non-reproducible run)"}\n')
        out.write(f'Clonotype Name Column,{clonotype_id_column}\n')
        out.write(f'Clonotype Names: Suffixed To Stay Unique,{len(name_collisions)}\n')
        out.write(f'Input TCRs: Unique After Deduplication,{tcr_list.shape[0]}\n')
        out.write(f'Input TCRs: Duplicates Removed,{removed_num}\n')
        out.write(f'Oligos: Failed To Generate,{invalid_cdr3_num}\n')

        out.write(f'Num Oligos: Pool A Under 350bp,{poolA_under350.shape[0]}\n')
        out.write(f'Num Oligos: Pool B Under 350bp,{poolB_under350.shape[0]}\n')
        out.write(f'Num Oligos: Pool A Over 350bp,{poolA_over350.shape[0]}\n')
        out.write(f'Num Oligos: Pool B Over 350bp,{poolB_over350.shape[0]}\n')

        out.write(f'Avg Oligo Size: Pool A Under 350bp,{poolA_under350["Length"].mean()}\n')
        out.write(f'Avg Oligo Size: Pool B Under 350bp,{poolB_under350["Length"].mean()}\n')
        out.write(f'Avg Oligo Size: Pool A Over 350bp,{poolA_over350["Length"].mean()}\n')        
        out.write(f'Avg Oligo Size: Pool B Over 350bp,{poolB_over350["Length"].mean()}\n')

def run_CLI():
    parser = argparse.ArgumentParser(description='TCRAFT-generate: Generate CDR3 Oligos from TCR sequence annotation data.')
    parser.add_argument('input_csv_path', type=str, help='Input CSV file with TCR information. Must contain the following columns: V_alpha, V_beta, CDR3_alpha, CDR3_beta, J_alpha, J_beta')
    parser.add_argument('--output_dir', type=str, help='Directory to save output files. Default is ./TCRAFT-generate_<date>', default=None)
    parser.add_argument('--clonotype_id_column', type=str, default=None, help=f'Input CSV column containing clonotype names. Auto-detected by default among: {", ".join(CLONOTYPE_CANDIDATE_COLUMNS)}. If the CSV holds none, names TCR_0001... are generated.')
    parser.add_argument('--seed', type=int, default=DEFAULT_SEED, help=f'Random seed for codon sampling. The same seed and input CSV always yield the same oligos. Default is {DEFAULT_SEED}.')
    parser.add_argument('--no_seed', action='store_true', help='Do not fix the random seed (non-reproducible run).')
    args = parser.parse_args()
    generate_cdr3_oligos(args.input_csv_path, args.output_dir, args.clonotype_id_column,
                         ref_dir.joinpath('random_1000nt.txt'),
                         seed=None if args.no_seed else args.seed)
    
if __name__ == "__main__":
    run_CLI()