from os import path,getcwd
from Bio import SeqIO,Entrez
import primer3
from pandas import DataFrame,__version__
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import sys
import numpy as np
import streamlit as st


@st.cache
def getpath():
    path_re = getcwd()
    return path_re

@st.cache
# Prepare system report
def get_system_report():
    report = {}
    report['primer_picker_version'] = "v0.0.1"
    report['python_version'] = sys.version[:5]
    report['pandas_version'] = __version__
    report['numpy_version'] = np.version.version
    return report

@st.cache
def seq_clean(raw):  # 不适合多个序列的情况
    raw = raw.upper()
    nucl = list('ATGC')
    out_seq = ''
    for nuc in raw:
        if nuc in nucl:
            out_seq += nuc
    return out_seq

@st.cache
def readFASTA(FASTA):
    Index = []
    Seq = []
    n = 0
    txt = ''
    for i in FASTA:
        if '>' in i and n == 0:
            txt = ''
            Index.append(i.replace('>', '').replace('\n', ''))  # replace('\n','')
            n += 1
        elif '>' in i and n > 0:
            Seq.append(txt)
            txt = ''
            Index.append(i.replace('>', '').replace('\n', ''))
        else:
            txt += i.strip()
    Seq.append(txt)
    return Seq, Index


def get_target(email: str, gene_id: dict, gene_name=''):
    Entrez.email = email
    if not gene_name:
        filename = gene_id['id']
    else:
        filename = gene_name
    if not path.isfile(filename):
        net_handle = Entrez.efetch(db=gene_id['db'], id=gene_id['id'], rettype=gene_id['rettype'],
                                   retmode=gene_id['retmode'])
        out_handle = open(filename, "w")
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print(filename + ' Data Saved.')
    else:
        print('Data found in given local address.')
    print('Parsing...')
    record = SeqIO.read(filename, 'gb')
    seq = str(record.seq)
    return record, seq


def design_primer(seq_args, global_args, out_address):
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args)
    info = '''Succeeded.  \n A total of {} primer pair(s) are designed.  \n For forward primer: {}  \n For reverse 
    primer: {}  \n For primer pairs: {}'''.format(
        primer3_result['PRIMER_PAIR_NUM_RETURNED'], primer3_result['PRIMER_LEFT_EXPLAIN'],
        primer3_result['PRIMER_RIGHT_EXPLAIN'], primer3_result['PRIMER_PAIR_EXPLAIN'])
    primer3_result_table_dict = {}
    for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        primer_id = str(i)
        for key in primer3_result:
            if '_' + primer_id + '_' in key:
                info_tag = key.replace("_" + primer_id, "")
                try:
                    primer3_result_table_dict[info_tag]
                except:
                    primer3_result_table_dict[info_tag] = []
                finally:
                    primer3_result_table_dict[info_tag].append(primer3_result[key])
    index = []
    for i in range(primer3_result["PRIMER_PAIR_NUM_RETURNED"]):
        index.append("PRIMER_PAIR_" + str(i + 1))
    primer3_result_df = DataFrame(primer3_result_table_dict, index=index)
    primer3_result_df = primer3_result_df.T
    primer3_result_df.to_csv(out_address)
    # print('Primer pair(s) Designed: \n', primer3_result_df)
    return primer3_result_df, info


def extract_pairs(primer3_result_df):
    # primer3_result_df = read_csv(out_address, header=0, index_col=0)
    left_primer = list(primer3_result_df.loc['PRIMER_LEFT_SEQUENCE'].values)
    right_primer = list(primer3_result_df.loc['PRIMER_RIGHT_SEQUENCE'].values)
    primer_list = []
    while left_primer:
        f = left_primer.pop(0)
        r = right_primer.pop(0)
        primer_list.append([f, r])
    return primer_list


def blastn(query_address: str, db_address: str, out_address1: str = '', evalue=0.001, identity=18, task='blastn',
           dust='yes'):
    blastn_cline = NcbiblastnCommandline(query=query_address, db=db_address, evalue=evalue, outfmt=5, out=out_address1,
                                         task=task, dust=dust)
    stout, stderr = blastn_cline()
    result_handle = open(out_address1)
    blast_record = NCBIXML.read(result_handle)
    e_value_thresh = evalue  # set E_value or other parameter and judge if exist
    identities = identity  # set identity for alignments,for primer design:length of primer-2 is recommended
    count = 0  # count number of blast hits
    name_list = []
    message = ''
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect <= e_value_thresh and hsp.identities >= identities:
                count += 1
                name_list.append(alignment.title)
                # 解决｜的问题
                sign = ''
                for a in hsp.match[0:75]:
                    if a == '|':
                        sign += '|'
                    else:
                        sign += ' '
                message += (
                    '  \n >> **sequence**: {0}  \n >> **length**: {1}  \n >> **identity**: {2}  \n >> **e-value**: {3}  '
                    '\n >>  `{4}`  '
                    '\n >>  `{5}`  '
                    '\n >>  `{6} `  \n >> '.format(
                        alignment.title, alignment.length, hsp.identities, hsp.expect, hsp.query[0:75] + '...',
                                                                                       sign + '...',
                                                                                       hsp.sbjct[0:75] + '...'))
    message += ' **{} similar sequence found.**  \n'.format(str(count))
    # st.write(message)
    return count, name_list, message

@st.cache
def set_args(job_name, target_seq, INCLUDED_REGION, state_target_region, TARGET_REGION, state_exclude_region,
             EXCLUDED_REGION_1, PRIMER_OPT_SIZE, PRIMER_PICK_INTERNAL_OLIGO, PRIMER_INTERNAL_MAX_SELF_END, primer_size,
             PRIMER_OPT_TM, tm, PRIMER_GC, PRIMER_MAX_POLY_X, PRIMER_INTERNAL_MAX_POLY_X, PRIMER_SALT_MONOVALENT,
             PRIMER_DNA_CONC, PRIMER_MAX_NS_ACCEPTED, PRIMER_MAX_SELF_ANY, PRIMER_MAX_SELF_END,
             PRIMER_PAIR_MAX_COMPL_ANY, PRIMER_PAIR_MAX_COMPL_END, PRIMER_PRODUCT_SIZE_RANGE, PRIMER_NUM_RETURN):
    seq_args = {  # Set parameters: name,seq, include region etc.
        'SEQUENCE_ID': job_name,
        'SEQUENCE_TEMPLATE': target_seq,
        'SEQUENCE_INCLUDED_REGION': [INCLUDED_REGION[0] - 1, INCLUDED_REGION[1] - INCLUDED_REGION[0]]
    }
    if state_target_region is True:
        seq_args['SEQUENCE_TARGET'] = [TARGET_REGION[0] - 1, TARGET_REGION[1] - TARGET_REGION[0]]
    else:
        seq_args['SEQUENCE_TARGET'] = []
    if state_exclude_region is True:
        seq_args['SEQUENCE_EXCLUDED_REGION'] = [[EXCLUDED_REGION_1[0] - 1, EXCLUDED_REGION_1[1] - EXCLUDED_REGION_1[0]]]
    else:
        seq_args['SEQUENCE_EXCLUDED_REGION'] = []

    global_args = {
        'PRIMER_OPT_SIZE': PRIMER_OPT_SIZE,
        'PRIMER_PICK_INTERNAL_OLIGO': PRIMER_PICK_INTERNAL_OLIGO,
        'PRIMER_INTERNAL_MAX_SELF_END': PRIMER_INTERNAL_MAX_SELF_END,
        'PRIMER_MIN_SIZE': primer_size[0],
        'PRIMER_MAX_SIZE': primer_size[1],
        'PRIMER_OPT_TM': PRIMER_OPT_TM,
        'PRIMER_MIN_TM': tm[0],
        'PRIMER_MAX_TM': tm[1],
        'PRIMER_MIN_GC': PRIMER_GC[0],
        'PRIMER_MAX_GC': PRIMER_GC[1],
        'PRIMER_MAX_POLY_X': PRIMER_MAX_POLY_X,
        'PRIMER_INTERNAL_MAX_POLY_X': PRIMER_INTERNAL_MAX_POLY_X,
        'PRIMER_SALT_MONOVALENT': PRIMER_SALT_MONOVALENT,
        'PRIMER_DNA_CONC': PRIMER_DNA_CONC,
        'PRIMER_MAX_NS_ACCEPTED': PRIMER_MAX_NS_ACCEPTED,
        'PRIMER_MAX_SELF_ANY': PRIMER_MAX_SELF_ANY,
        'PRIMER_MAX_SELF_END': PRIMER_MAX_SELF_END,
        'PRIMER_PAIR_MAX_COMPL_ANY': PRIMER_PAIR_MAX_COMPL_ANY,
        'PRIMER_PAIR_MAX_COMPL_END': PRIMER_PAIR_MAX_COMPL_END,
        'PRIMER_PRODUCT_SIZE_RANGE': [PRIMER_PRODUCT_SIZE_RANGE],
        'PRIMER_NUM_RETURN': PRIMER_NUM_RETURN
    }
    return seq_args, global_args


