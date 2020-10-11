##########################
# written by Hasan RH Al-Saedi and Sébastien Leblanc 
# under the supervision of Xavier Roucou and Marie A Brunet. The code is used in the OpenProt pipeline, see: Brunet MA et al, 2021
##########################

import os
import sys
import re
import argparse
import itertools as itt
import numpy as np
from multiprocessing import Pool

parser = argparse.ArgumentParser(description='Compute theoretical peptides and detected coverage.')
parser.add_argument('file_path', type=str, help='Path to text file containing the list of protein sequences.')
parser.add_argument('proteome_tsv', type=str, help='Path to text file containing the full proteome for unicity checking.')
parser.add_argument('out_dir', type=str, help='Directory for output of tsv.')
parser.add_argument('--n_cpu', type=int, help='Number of cpu for paralell processing.', default=4)
parser.add_argument('--progress_iter_size', type=int, help='Number of lines to report progress at.', default=10000)
args = parser.parse_args()

aa_weights = {
    'X': 110, 'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'E': 147, 'Q': 146,
    'G': 75, 'H': 155, 'I': 131, 'L':131,'K': 146, 'M': 149, 'O' : 255, 'F': 165, 'P': 115,
    'S': 105, 'T': 119, 'W': 204, 'Y': 181,  'U':168, 'V': 117, '1': 131,
}
columns = ['prot_seq_id', 'gene_symbol', 'tt_type', 'prot_seq', 'pep_seq_ls']

def load_tsv(file_path):
    ''' IN: tab delimited text file with columns as defined here
        OUT: list of dicts with columns as keys and text as values

    '''
    prot_dicts = []
    with open(file_path, 'r') as f:
        for l in f:
            ls = l.strip().split('\t')
            line = dict(zip(columns, ls))
            if any(x not in aa_weights for x in line['prot_seq']): continue # skip proteins with unkown aa in the sequence i.e. "*"
            if 'pep_seq_ls' not in line: # no peptides detected, append empty list
                line['pep_seq_ls'] = []
            else:
                line['pep_seq_ls'] = line['pep_seq_ls'].split(',')
            prot_dicts.append(line)
    return prot_dicts

def trypsin_digest(protseq, debug=False):
    ### mc -> miss cleavage
    raw_peps = []
    raw_peps_0mc = re.sub("K(?!P)", "K,", protseq)
    raw_peps_0mc = re.sub("R(?!P)", "R,", raw_peps_0mc)
    raw_peps_0mc = raw_peps_0mc.split(",")
    raw_peps_1mc = [raw_peps_0mc[n] + raw_peps_0mc[n+1] for n in range(len(raw_peps_0mc) - 1)]
    raw_peps_2mc = [raw_peps_1mc[n] + raw_peps_0mc[n+2] for n in range(len(raw_peps_1mc) - 1)]
    raw_peps.extend(raw_peps_0mc)
    raw_peps.extend(raw_peps_1mc)
    raw_peps.extend(raw_peps_2mc)
    #if debug: log('{} >> {} >> {}'.format(raw_peps_0mc, raw_peps_1mc, raw_peps_2mc))
    return raw_peps

def pep_size_filter(raw_peps):
    pred_peps = [raw_pep for raw_pep in raw_peps if len(raw_pep) >= 7 and sum(aa_weights[aa] for aa in raw_pep) < 4600]
    return pred_peps

def unify_isobaric_aas(aa_seq):
    common_char = '1'
    return aa_seq.replace('I', common_char).replace('L', common_char)

def check_pep_unicity(pep, prot_dict):
    if prot_dict['tt_type'] == 'reference': # peptide unicity rules for refs (only compared to refs of other genes..)
        check_seqs = [p['prot_seq'] for p in full_proteome if p['gene_symbol'] != prot_dict['gene_symbol'] and p['tt_type']=='reference']
    else: # peptide unicity rules for alts (and iso) *** verify
        check_seqs = [p['prot_seq'] for p in full_proteome if p['tt_type'] == 'reference']
        check_seqs += [p['prot_seq'] for p in full_proteome if p['gene_symbol'] != prot_dict['gene_symbol'] and any(x in p['tt_type'] for x in ['cryptic', 'isoform']) ]

    pep = unify_isobaric_aas(pep)
    for protseq in list(set(check_seqs)):
        if pep in unify_isobaric_aas(protseq): return False

    return True

def map_peptides(protseq, peps):
    peps_locs = []
    for pep in peps:
        if pep not in protseq: continue
        loc_st  = protseq.find(pep)
        loc_end = loc_st + len(pep)
        pep_loc = protseq[loc_st:loc_end]
        peps_locs.append((pep_loc, loc_st, loc_end))
    return peps_locs

def compute_pep_coverage(protseq, peps_locs):
    protvec = [0] * len(protseq)
    for pepseq_start_end in peps_locs:
        pepseq, start, end = pepseq_start_end
        protvec[start:end] = [1] * len(pepseq)
    coverage = (sum(protvec) / len(protseq))*100
    return coverage

def get_theo_peps(prot_dict):
    tryp_peps = pep_size_filter(trypsin_digest(prot_dict['prot_seq']))
    unique_tryp_peps = [pep for pep in tryp_peps if check_pep_unicity(pep, prot_dict)]
    return unique_tryp_peps

def compute_theo_coverage(prot_dict):
    unique_tryp_peps = get_theo_peps(prot_dict)
    peps_locs = map_peptides(prot_dict['prot_seq'], unique_tryp_peps)
    coverage = compute_pep_coverage(prot_dict['prot_seq'], peps_locs)
    return coverage

def compute_detect_coverage(prot_dict):
    peps_locs = map_peptides(prot_dict['prot_seq'], prot_dict['pep_seq_ls'])
    coverage = compute_pep_coverage(prot_dict['prot_seq'], peps_locs)
    return coverage

def compute_all_coverages(prot_dict):
    prot_dict['theo_coverage'] = compute_theo_coverage(prot_dict)
    prot_dict['detect_coverage'] = compute_detect_coverage(prot_dict)
    return prot_dict

def dict_list_to_tsv(dict_ls, file_path):
    write_cols = columns + ['theo_coverage', 'detect_coverage']
    with open(file_path, 'w') as f:
        for line in dict_ls:
            line['pep_seq_ls'] = ','.join(line['pep_seq_ls'])
            l = '\t'.join([str(line[c]) for c in write_cols]) + '\n'
            f.write(l)

if __name__ == '__main__':
    args = vars(args)

    prot_dicts = load_tsv(args['file_path'])
    print('Processing {} proteins.'.format(len(prot_dicts)))

    full_proteome = load_tsv(args['proteome_tsv'])

    pool = Pool(processes=args['n_cpu'])
    progress = 0
    results = []
    for i, r in enumerate(pool.imap_unordered(compute_all_coverages, prot_dicts)):
        if progress == args['progress_iter_size']:
            progress = 0
            print('completion: {}% ({} prots)'.format((i/len(prot_dicts))*100, i))
        progress += 1
        results.append(r)
    pool.close()
    pool.join()

    output_file = args['out_dir'] + '/' + os.path.basename(args['file_path'])  +'.out.tsv'
    dict_list_to_tsv(results, output_file)
    print('Done. Outputed to {}'.format(output_file))