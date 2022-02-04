import streamlit as st
from Basic_Functions.Primer_picker_funtion import *
# from dna_features_viewer import GraphicFeature, GraphicRecord
from os import path,remove
# from Bio.Seq import Seq
# import matplotlib.pyplot as plt



def primer_picker_single_main():
    with st.sidebar:
        # Define sidebar of primer picker
        st.header('Primer Picker settings')
        db = st.selectbox('Please Choose database to blast from',
                                  ['B. napus', 'B. rapa', 'B. oleracea', 'R. sativus'])
        if db not in ['B. napus', 'R. sativus']:
            st.warning('Currently not available.')

        out_address = r'./primer_picker_cache/Primer_Picker_Output/'
        evalue = st.number_input('E-value:', value=1.0)
        mismatch = st.number_input('Mismatch for blast:', value=2)
        state_button_primer = st.button('Get Primers only')
        state_button_blast = st.button('Blast Primers')

    # main page info
    col1, col2 = st.columns(2)
    with col1:
        job_name = st.text_input('Job Name:', 'AHK2_A03')
    with col2:
        email = st.text_input('Please enter e-mail address:', 'zhu.xingyu.r1@dc.tohoku.ac.jp')

    sample_seq = '''>sample_seq\nggggtcttcccccaagaattggtagatgtggttttggtgccggccgacatgtctatagcttgtgagctctcgaatcccaatttaaagaagacaaaagcagaaaaacggataccaaccaagattctgcttatccgggttttatgtggcttagttgttctctggctctgcttaagcttaagcttaggcttcttatgtatatgcaagaagaaagaggctgctgctgctgctgatgactcttcttcttctgctaaagggatgttgttcaggaatcagagcagaagtgagattgatgctatgctttctctcttctttgattcaaatcaggtttcccactacttctttcttctttttattattgtttgacaatttaatttacttaaagcttttgaataattggagaacaggtaacatcttttgaatgtcgcaaggagaatggtggtattacatgttccttgtcaacacgttccgagaaaggagacgaggaggaggaggaggctaagagacatgttgttgcagaggttatgtcatcatctgagaatgaagaagaaggaggtgtgctgcatcaggttgtgttgttttatgtaatgaacaaatggcattggtggttggtcctttgtgtactactagtgggcggcggccgtgtgatatttgtaagaaaagaagtttcttctttagtacaagacaagcagcagcagcagcagcaatgcaaaacagctgggaagtggaggaagaacatgcttctactcggcatcatcgcgggagtttctttgtctgttttatggttttgggataccaacgagaagatcttgttccaaaggaaagagacgttaaccaacatgtgtgaggaacgagctcgggtgttgcaggaccagttcaatgttagcatgaaccatgtccacgccttgtccatcctcgtctctacctttcaccatggcaaaaccccttctgccattgatcaggtaaactaagagtccgaatggtaacatgtggtgcggttctaatagtaacaaaatcgtatagatacatgtgtaaatctagagattttttaatgttaattgcggtgcgattctaatagtaaaaaaaagtatagatacatgtgtatatctagagatttttgttattgttaactacgatggtgttttaatagtaacaaaaacgtatagatacatatgtatatttagagatttttgttaatgttaattgcggcgcggttctaatagtaaccaaaaagtatatatacatgtgtatatctagagatttttttaactgcggtgctgttccaatagtaacaaaaacgtatagatacatgtgtaaatctagagatttttgttaatgttaattgtggtgcagttctaatagtaacaaaaacgtatagatacatgtgtatatctaaagatttttgttaatgttaattgtggtgcggttccaatagtaacaaaaacgtatagatacatgtgtaaatctagagatttttgttaatgttaattgtggtgcggttctatagtaacaaaaacgtatagatacatgtgtatatctagagatttttgttaatgctaattgtgctgcggtgcggttagaaataaaattccattcttaagatttttaactgttgattttttttttttcag'''

    fasta = st.text_area('Input fasta sequence.', sample_seq, height=250)

    raw_seq, index = readFASTA(fasta)

    target_seq = seq_clean(raw_seq[0])

    state_target_region = st.checkbox('Target region needed.')
    state_exclude_region = st.checkbox('Exclude region needed.')
    # PCR product parameters
    a = st.container()
    a.header('PCR Product Parameters')

    col3, col4 = a.columns(2)
    with col3:
        start = int(
            st.number_input('PCR product size from:', min_value=1, max_value=3000, value=100, key='product_size_start'))
    with col4:
        end = int(st.number_input('to', min_value=1, max_value=3000, value=500, key='product_size_end'))
    PRIMER_PRODUCT_SIZE_RANGE = [start, end]

    INCLUDED_REGION = a.slider('Select Region to Include', min_value=1, max_value=len(target_seq) + 1,
                               value=(1, len(target_seq) + 1), step=10,
                               key='include_region')
    if state_target_region is True:
        TARGET_REGION = a.slider('Select Region as target', 1, len(target_seq) + 1, (1, len(target_seq) + 1), step=10,
                                 key='target_region')
    else:
        TARGET_REGION = False
    if state_exclude_region is True:
        EXCLUDED_REGION_1 = a.slider('Select Region to exclude', 1, len(target_seq) + 1, (1, len(target_seq) + 1),
                                     step=10,
                                     key='exclude_region_1')
    else:
        EXCLUDED_REGION_1 = False

    # Primer parameters
    c = st.container()
    c.header('Primer Parameters')
    col5, col6, col7 = c.columns([1, 1, 6])
    # with col5:
    #     st.markdown('&nbsp ')
    #     st.markdown('##### {}    &nbsp    &nbsp &nbsp '.format('Primer Size:'))
    #     st.markdown('##### {}    &nbsp    &nbsp &nbsp '.format('Tm Value:'))
    #     st.markdown('##### {}    &nbsp    &nbsp &nbsp '.format('Number of Primer:'))

    with col6:
        PRIMER_OPT_SIZE = int(st.number_input('Primer Optimal Size:', 20, key='PRIMER_OPT_SIZE'))
        PRIMER_OPT_TM = int(st.number_input('Primer Optimal Tm:', 60, key='PRIMER_OPT_TM'))
    with col7:
        primer_size = st.slider('Select Range of Primer Length', 15, 40, (18, 25),
                                key='target_region')

        tm = st.slider('Select Range of Tm Value', 50, 70, (58, 62),
                       key='tm_value')

    PRIMER_GC = c.slider('Select Range of GC content', 20, 80, (40, 60),
                         key='gc_content')

    PRIMER_NUM_RETURN = int(
        c.number_input('Number of primers to return:', min_value=1, max_value=200, value=20,
                       key='PRIMER_NUM_RETURN',
                       help='Number of primers to return'))

    # Advanced Primer Parameters
    with st.expander('Advanced Primer Parameters', False):
        state_ambiguouss = st.checkbox("Remove ambiguous primer. ", key='state_ambiguous')
        PRIMER_MAX_POLY_X = int(st.number_input('PRIMER_MAX_POLY_X:', 100, key='PRIMER_MAX_POLY_X'))

        internal_oligo = st.checkbox('Pick Internal Primer', key='internal_oligo')
        if internal_oligo:
            PRIMER_PICK_INTERNAL_OLIGO = 1
            PRIMER_INTERNAL_MAX_SELF_END = int(
                st.number_input('PRIMER_INTERNAL_MAX_SELF_END:', 8, key='internal_self_end'))
            PRIMER_INTERNAL_MAX_POLY_X = int(
                st.number_input('PRIMER_INTERNAL_MAX_POLY_X:', 100, key='PRIMER_INTERNAL_MAX_POLY_X'))
        else:
            PRIMER_PICK_INTERNAL_OLIGO = 0
            PRIMER_INTERNAL_MAX_SELF_END = 8
            PRIMER_INTERNAL_MAX_POLY_X = 100

        PRIMER_SALT_MONOVALENT = int(st.number_input('PRIMER_SALT_MONOVALENT:', 50, key='PRIMER_SALT_MONOVALENT'))
        PRIMER_DNA_CONC = int(st.number_input('PRIMER_DNA_CONC:', 50, key='PRIMER_DNA_CONC'))
        PRIMER_MAX_NS_ACCEPTED = int(st.number_input('PRIMER_MAX_NS_ACCEPTED:', 0, key='PRIMER_MAX_NS_ACCEPTED'))
        PRIMER_MAX_SELF_ANY = int(st.number_input('PRIMER_MAX_SELF_ANY:', 12, key='PRIMER_MAX_SELF_ANY'))
        PRIMER_MAX_SELF_END = int(st.number_input('PRIMER_MAX_SELF_END:', 8, key='PRIMER_MAX_SELF_END'))
        PRIMER_PAIR_MAX_COMPL_ANY = int(
            st.number_input('PRIMER_PAIR_MAX_COMPL_ANY:', 12, key='PRIMER_PAIR_MAX_COMPL_ANY'))
        PRIMER_PAIR_MAX_COMPL_END = int(
            st.number_input('PRIMER_PAIR_MAX_COMPL_END:', 8, key='PRIMER_PAIR_MAX_COMPL_END'))

    if state_button_primer is True:
        seq_args, global_args = set_args(job_name, target_seq, INCLUDED_REGION, state_target_region, TARGET_REGION,
                                         state_exclude_region, EXCLUDED_REGION_1, PRIMER_OPT_SIZE,
                                         PRIMER_PICK_INTERNAL_OLIGO, PRIMER_INTERNAL_MAX_SELF_END, primer_size,
                                         PRIMER_OPT_TM, tm, PRIMER_GC, PRIMER_MAX_POLY_X,
                                         PRIMER_INTERNAL_MAX_POLY_X, PRIMER_SALT_MONOVALENT, PRIMER_DNA_CONC,
                                         PRIMER_MAX_NS_ACCEPTED, PRIMER_MAX_SELF_ANY, PRIMER_MAX_SELF_END,
                                         PRIMER_PAIR_MAX_COMPL_ANY, PRIMER_PAIR_MAX_COMPL_END,
                                         PRIMER_PRODUCT_SIZE_RANGE, PRIMER_NUM_RETURN)
        primer_df, message = design_primer(seq_args, global_args,
                                           out_address + job_name + '_Original_Primer.txt')  # To get all primer pairs in the form of dataframe
        st.info(message)
        # show dataframe
        primer_df = primer_df.astype(str)
        st.dataframe(primer_df)

    if state_button_blast:
        db_dic = {
            'B. napus': r'./primer_picker_cache/db/B.napus/Brassica_napus_v4.1.chromosomes',
            'B. rapa': '',
            'B. oleracea': '',
            'R. sativus': r'./primer_picker_cache/db/R.sativus/RSAskr_r1.0'
        }
        db_address = db_dic[db]

        seq_args, global_args = set_args(job_name, target_seq, INCLUDED_REGION, state_target_region, TARGET_REGION,
                                         state_exclude_region,
                                         EXCLUDED_REGION_1, PRIMER_OPT_SIZE, PRIMER_PICK_INTERNAL_OLIGO,
                                         PRIMER_INTERNAL_MAX_SELF_END, primer_size,
                                         PRIMER_OPT_TM, tm, PRIMER_GC, PRIMER_MAX_POLY_X, PRIMER_INTERNAL_MAX_POLY_X,
                                         PRIMER_SALT_MONOVALENT,
                                         PRIMER_DNA_CONC, PRIMER_MAX_NS_ACCEPTED, PRIMER_MAX_SELF_ANY,
                                         PRIMER_MAX_SELF_END,
                                         PRIMER_PAIR_MAX_COMPL_ANY, PRIMER_PAIR_MAX_COMPL_END,
                                         PRIMER_PRODUCT_SIZE_RANGE, PRIMER_NUM_RETURN)
        st.header('Primer Pick Result')
        primer_df, message = design_primer(seq_args, global_args,
                                           out_address + job_name + '_Original_Primer.txt')  # To get all primer pairs in the form of dataframe
        st.info(message)
        primer_list = extract_pairs(primer_df)  # Turn dataframe into list of primer pairs
        # show dataframe
        primer_df = primer_df.astype(str)
        st.dataframe(primer_df)

        # st.warning('Primer designed, checking for specificity...')
        list_result = []
        count = 0
        current = 0
        list_r = []
        list_f = []
        message = ''

        my_bar = st.progress(0)
        total = len(primer_list)

        if state_ambiguouss is False:
            for pair in primer_list:
                current += 1
                percent = current / total
                my_bar.progress(percent)
                short = ''
                short += '  \n  \n #### Checking Primer pair {} :  \n '.format(str(count + 1))
                forward = pair[0]
                tmp_f = open(out_address + "{}_f.fasta".format(job_name), "w")
                tmp_f.write("%s" % forward)
                tmp_f.close()
                reverse = pair[1]
                tmp_r = open(out_address + "{}_r.fasta".format(job_name), "w")
                tmp_r.write("%s" % reverse)
                tmp_r.close()
                short += '- **Left Primer**:  \n '
                f, f_n, f_short = blastn(out_address + "{}_f.fasta".format(job_name), db_address=db_address,
                                         out_address1=out_address + job_name + '_PBL.xml',
                                         evalue=evalue, identity=len(forward) - mismatch, task='blastn-short',
                                         dust='no')
                short += f_short
                short += '- **Right Primer**:  \n '
                r, r_n, r_short = blastn(out_address + "{}_r.fasta".format(job_name), db_address=db_address,
                                         out_address1=out_address + job_name + '_PBL.xml',
                                         evalue=evalue, identity=len(reverse) - mismatch, task='blastn-short',
                                         dust='no')
                short += r_short
                message += short + '  \n *** \n '
                list_f.append(f)
                list_r.append(r)

                if r == 1 or f == 1:
                    list_result.append(count)
                    # st.write('Primer Pair '+str(count)+': Forward/Reverse primer found '+str(f)+'/'+str(r)+' products on intended
                    # targets,respectively.')
                count += 1

        else:
            for pair in primer_list:
                current += 1
                percent = current / total
                my_bar.progress(percent)
                short = ''
                short += '  \n  \n #### Checking Primer pair {} :  \n '.format(str(count + 1))
                forward = pair[0]
                tmp_f = open(out_address + "{}_f.fasta".format(job_name), "w")
                tmp_f.write("%s" % forward)
                tmp_f.close()
                reverse = pair[1]
                tmp_r = open(out_address + "{}_r.fasta".format(job_name), "w")
                tmp_r.write("%s" % reverse)
                tmp_r.close()
                short += '- **Left Primer**:  \n '
                f, f_n, f_short = blastn(out_address + "{}_f.fasta".format(job_name), db_address=db_address,
                                         out_address1=out_address + job_name + '_PBL.xml',
                                         evalue=evalue, identity=len(forward) - mismatch, task='blastn-short',
                                         dust='no')
                short += f_short
                short += '- **Right Primer**:  \n '
                r, r_n, r_short = blastn(out_address + "{}_r.fasta".format(job_name), db_address=db_address,
                                         out_address1=out_address + job_name + '_PBL.xml',
                                         evalue=evalue, identity=len(reverse) - mismatch, task='blastn-short',
                                         dust='no')
                short += r_short
                message += short + '  \n *** \n '
                list_f.append(f)
                list_r.append(r)

                if r == 1 and f == 1:
                    list_result.append(count)
                    # st.write('Primer Pair '+str(count)+': Forward/Reverse primer found '+str(f)+'/'+str(r)+' products on intended
                    # targets,respectively.')
                count += 1

        with st.expander('Show Blast Log.'):
            st.markdown(message)

        primer_df.loc['LEFT_PRIMER_PRODUCTS'] = list_f
        primer_df.loc['RIGHT_PRIMER_PRODUCTS'] = list_r
        pd_result = primer_df.iloc[:, list_result]
        if pd_result.empty:
            st.error('\nNo specified primer pair found. Try change parameters or increase number of primers to return.')
        else:
            st.subheader('\nResult of useful primer(s):\n')
            pd_result = pd_result.astype(str)
            st.dataframe(pd_result)
            st.balloons()
            done = True
            st.success('Blast succeeded.')

        if path.exists(out_address + "{}_f.fasta".format(job_name)):
            remove(out_address + "{}_f.fasta".format(job_name))
        if path.exists(out_address + "{}_r.fasta".format(job_name)):
            remove(out_address + "{}_r.fasta".format(job_name))
        if path.exists(out_address + job_name + '_PBL.xml'):
            remove(out_address + job_name + '_PBL.xml')