from Basic_Functions.Primer_picker_funtion import *


def s_haplotype_main():
    # sidebar
    with st.sidebar:
        # Define sidebar of primer picker
        st.header('S-haplotype blast settings')
        blast_type = st.radio('Please Choose type of blast to perform:',
                              ['Nucleotide Blast', 'Protein Blast', 'Nucl -> Protein', 'Protein -> Nucl'])
        if blast_type not in ['Nucleotide Blast', 'Protein Blast']:
            st.warning('Currently not available.')
        db = st.multiselect('Database to blast from:',
                            options=['B.rapa', 'B.oleracea', 'R.sativus', 'R.raphanistrum'],
                            default=['B.rapa', 'B.oleracea', 'R.sativus', 'R.raphanistrum'])
        st.write(db,type(db))
        db_dic = {
            'Nucleotide Blast': r'./primer_picker_cache/db/B.napus/Brassica_napus_v4.1.chromosomes',
            'Protein Blast': '',
            'Nucl -> Protein': '',
            'Protein -> Nucl': ''
        }

        db_address = db_dic[db]
        out_address = r'./primer_picker_cache/Primer_Picker_Output/'
        evalue = st.number_input('E-value:', value=1.0)
        mismatch = st.number_input('Mismatch for blast:', value=2)
        state_button_blast = st.button('Blast')

    # main function
    job_name = st.text_input('Job Name:', value='Default')
    fasta = st.text_area('Input fasta sequence.',
                         value='>test_seq\nTGAACAGATTCCTTACATCATCGAGAAATTTCGACGATCCATCAAGCGGGGATTACTCGTACAA',
                         height=250)
    raw_seq, index = readFASTA(fasta)

    if state_button_blast is True:
        st.balloons()
        pass
