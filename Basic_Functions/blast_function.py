import streamlit as st
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML


def show_align(query, sbjct, sign):
    out = ''
    count = 1
    while len(query) > 120:
        out += ' \n >>  `query:  {0}  {3}~{4}`   \n  >>  `match:  {1}  {3}~{4}`  \n  >>  `sbjct:  {2}  {3}~{4}` \n'.format(
            query[0:120], sign[0:120], sbjct[0:120], (count - 1) * 120 + 1, count * 120)
        query = query[120:]
        sbjct = sbjct[120:]
        sign = sign[120:]
        count += 1
    else:
        out += ' \n >>  `query:  {0}  {3}~{4}`   \n  >>  `match:  {1}  {3}~{4}`  \n  >>  `sbjct:  {2}  {3}~{4}` \n'.format(
            query, sign, sbjct, (count - 1) * 120 + 1, (count - 1) * 120 + len(query))
    return out


def get_db_address(blast_type, ls):
    address = r'./primer_picker_cache/db/s_haplotype/full_db/'

    if len(ls) == 4:
        address += 'all_'
    else:
        if 'B.rapa' in ls:
            address += 'Br_'
        if 'B.oleracea' in ls:
            address += 'Bo_'
        if 'R.sativus' in ls:
            address += 'Rs_'
        if 'R.raphanistrum' in ls:
            address += 'Rr_'
    if blast_type in ['Nucleotide Blast', 'Protein -> Nucl']:
        address += 'nuc'
    else:
        address += 'pep'
    return address


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
                sign = ''
                for a in hsp.match:
                    if a == '|':
                        sign += '|'
                    else:
                        sign += ' '
                out_align = show_align(hsp.query, hsp.sbjct, sign)
                message += (
                    '  \n >> **sequence**: {0}  \n >> **length**: {1}  \n >> **identity**: {2} out of {3}  \n >> **e-value**: {4}  '
                    '{5} \n >> '.format(
                        alignment.title, alignment.length, hsp.identities, hsp.align_length, hsp.expect, out_align))
    message += ' ⇨ **{} similar sequence found.**  \n'.format(str(count))
    # st.write(message)
    return count, name_list, message


def blastp(query_address: str, db_address: str, out_address1: str = '', evalue=0.001, identity=18, task='blastp'):
    blastp_cline = NcbiblastpCommandline(query=query_address, db=db_address, evalue=evalue, outfmt=5, out=out_address1,
                                         task=task)
    stout, stderr = blastp_cline()
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
                sign = ''
                for a in hsp.match:
                    if a == ' ':
                        sign += '-'
                    else:
                        sign += '|'
                out_align = show_align(hsp.query, hsp.sbjct, sign)
                message += (
                    '  \n >> **sequence**: {0}  \n >> **length**: {1}  \n >> **identity**: {2} out of {3}  \n >> **e-value**: {4}  '
                    '{5} \n >> '.format(
                        alignment.title, alignment.length, hsp.identities, hsp.align_length, hsp.expect, out_align))
    message += ' ⇨ **{} similar sequence found.**  \n'.format(str(count))
    # st.write(message)
    return count, name_list, message


def blastx(query_address: str, db_address: str, out_address1: str = '', evalue=0.001, identity=18, task='blastx'):
    blastx_cline = NcbiblastxCommandline(query=query_address, db=db_address, evalue=evalue, outfmt=5, out=out_address1,
                                         task=task)
    stout, stderr = blastx_cline()
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
                sign = ''
                for a in hsp.match:
                    if a == ' ':
                        sign += '-'
                    else:
                        sign += '|'
                out_align = show_align(hsp.query, hsp.sbjct, sign)
                message += (
                    '  \n >> **sequence**: {0}  \n >> **length**: {1}  \n >> **identity**: {2} out of {3}  \n >> **e-value**: {4}  '
                    '{5} \n >> '.format(
                        alignment.title, alignment.length, hsp.identities, hsp.align_length, hsp.expect, out_align))
    message += ' ⇨ **{} similar sequence found.**  \n'.format(str(count))
    # st.write(message)
    return count, name_list, message


def tblastn(query_address: str, db_address: str, out_address1: str = '', evalue=0.001, identity=18, task='tblastn'):
    tblastn_cline = NcbitblastnCommandline(query=query_address, db=db_address, evalue=evalue, outfmt=5,
                                           out=out_address1,
                                           task=task)
    stout, stderr = tblastn_cline()
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
                sign = ''
                for a in hsp.match:
                    if a == ' ':
                        sign += '-'
                    else:
                        sign += '|'
                out_align = show_align(hsp.query, hsp.sbjct, sign)
                message += (
                    '  \n >> **sequence**: {0}  \n >> **length**: {1}  \n >> **identity**: {2} out of {3}  \n >> **e-value**: {4}  '
                    '{5} \n >> '.format(
                        alignment.title, alignment.length, hsp.identities, hsp.align_length, hsp.expect, out_align))
    message += ' ⇨ **{} similar sequence found.**  \n'.format(str(count))
    # st.write(message)
    return count, name_list, message
