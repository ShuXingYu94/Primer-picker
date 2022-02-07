from Basic_Functions.primer_picker_funtion import *
from Basic_Functions.blast_function import *
from os import path, remove


@st.cache
def get_sample_seq(type_blast):
    if type_blast in ['Nucleotide Blast', 'Nucl -> Protein']:
        out = '>test_seq\nATGCAAGGTGTACGATACATCTATCACCATTTTTACACCTCCTTGCTCGTCTTCGTTGTCATGATTCTATTTCGTTCTGCCCTTTCGA' \
              'TCTATATCAACACTTTGTCCTCTACAGAATCTCTTACAATCTCAAACAACAGAACACTTGTATCTCCCGGTGATGTTTTCGAGCTCGGTTTCTTCACTC' \
              'CCGGATCAAGTTCTCGTTGGTATCTCGGGATATGGTACAAGAAACTGCCCTATATAACCTATGTATGGGTTGCCAAC\n>test_seq_2\nCAAGATC' \
              'TTTATGTCAGATTGGCTGCCGCTGATCTTGTTAAAAAGAGAGACGCGAATTGGAAAATCATAAGTTTGATTGTTGGAGTTAGTGTTGTTCTGCTTCTTA' \
              'TGATCATGTTCTGCCTTTGGAAAAAGAAACAAAATCGAGCAAAAGCAATGGCATCATCTATTGTAAATCACCAGAGAAACCAAAATGTACTTATGAACA' \
              'CGATGACACAATCAAACAAGAGACAGTTGTCTAGAGAGAACAAAATTGAGGAAT'
    else:
        out = '>test_pep\nMQGVRYIYHHFYTSLLVFVVMILFRSALSIYINTLSSTESLTISNNRTLVSPGDVFELGFFTPGSSSRWYLGIWYKKLPYITYVWVANR' \
              'DNPLSNSTGTLKISGNNLFLLGDSNKSIWSTNLTRGNERSPVVAELLANGNFVMRDSNNNDASGFLWQSFDYPTDTLLPEMKLGYDLKTGLNRFLTSSRN' \
              'FDDPSSGDYSYKLEPRRLPEFYLLLGDVREHRSGPWNGIQFSGIPEDQKLSYMVYNFTKNSEEVAYTFRMTNNSFYSRLTINSE\n>test_pep_2\nQD' \
              'LYVRLAAADLVKKRDANWKIISLIVGVSVVLLLMIMFCLWKKKQNRAKAMASSIVNHQRNQNVLMNTMTQSNKRQLSRENKIEEFELPLIELEAVVKATE' \
              'NFSNCNELGRSGFGIVYKGMLDGQEVAVKRLSKTSLQGIDEFMNEVRLIARLQHINLVRILGCCIEADEKILIYEYLENSSLDYFLFGKKRSSNLNWKDR' \
              'FAITNGVARGLLYLHQDSRFRIIHRDLKPGNILLDKYMIPKISDFGMARIFARDETQARTDNAVGTYGYMSPEY'
    return out


def s_haplotype_main():
    # sidebar
    with st.sidebar:
        # Define sidebar of primer picker
        st.header('S-haplotype blast settings')
        blast_type = st.radio('Please Choose type of blast to perform:',
                              ['Nucleotide Blast', 'Protein Blast', 'Nucl -> Protein', 'Protein -> Nucl'])
        # if blast_type not in ['Nucleotide Blast', 'Protein Blast']:
        #     st.warning('Currently not available.')
        db = st.multiselect('Database to blast from:',
                            options=['B.rapa', 'B.oleracea', 'R.sativus', 'R.raphanistrum'],
                            default=['B.rapa', 'B.oleracea', 'R.sativus', 'R.raphanistrum'])
        out_address = r'./primer_picker_cache/S_haplotype_output/'
        state_button_blast = st.button('Blast')

    # main function
    sample_seq = get_sample_seq(blast_type)
    fasta = st.text_area('Input fasta sequence.',
                         value=sample_seq,
                         height=230, key='seq')
    raw_seq, index = readFASTA(fasta)
    job_name = st.text_input('Job Name:', value=str(index[0]))
    evalue = st.number_input('E-value:', value=0.000001, step=0.000001, key='evalue', format="%.6f")
    identity = st.number_input('Percentage of identity (%):', value=95, key='identity')
    mismatch = st.number_input('Number of mismatch:', value=20, key='mismatch')
    with st.spinner('Please wait until blast complete...'):
        if state_button_blast is True:
            db_address = get_db_address(blast_type, db)
            message = ''

            # Progress bar visualize
            current = 0
            my_bar = st.progress(0)
            total = len(raw_seq)
            for n in range(0, len(raw_seq)):
                current += 1
                percent = current / total
                my_bar.progress(percent)

                seq = raw_seq[n]
                id_1 = len(seq) * identity / 100
                id_2 = len(seq) - mismatch
                if id_1 > id_2:
                    iden = id_1
                else:
                    iden = id_2
                # st.write(id_1)
                # st.write(id_2)

                short = ''
                short += '  \n  \n #### Blast Result for {} :  \n '.format(index[n])
                tmp_seq = open(out_address + "{}.fasta".format(job_name), "w")
                tmp_seq.write("%s" % seq)
                tmp_seq.close()
                # short += '- **Left Primer**:  \n '
                if blast_type == 'Nucleotide Blast':
                    if len(seq) < 50:
                        task = 'blastn-short'
                    else:
                        task = 'blastn'
                    target, target_n, target_short = blastn(out_address + "{}.fasta".format(job_name),
                                                            db_address=db_address,
                                                            out_address1=out_address + job_name + '_SBL.xml',
                                                            evalue=evalue, identity=iden,
                                                            task=task,
                                                            dust='no')
                elif blast_type == 'Protein Blast':
                    if len(seq) < 30:
                        task = 'blastp-short'
                    else:
                        task = 'blastp'
                    target, target_n, target_short = blastp(out_address + "{}.fasta".format(job_name),
                                                            db_address=db_address,
                                                            out_address1=out_address + job_name + '_SBL.xml',
                                                            evalue=evalue, identity=iden,
                                                            task=task)
                elif blast_type == 'Nucl -> Protein':
                    target, target_n, target_short = blastx(out_address + "{}.fasta".format(job_name),
                                                            db_address=db_address,
                                                            out_address1=out_address + job_name + '_SBL.xml',
                                                            evalue=evalue, identity=iden / 3,
                                                            task='blastx')
                elif blast_type == 'Protein -> Nucl':
                    target, target_n, target_short = tblastn(out_address + "{}.fasta".format(job_name),
                                                             db_address=db_address,
                                                             out_address1=out_address + job_name + '_SBL.xml',
                                                             evalue=evalue, identity=iden,
                                                             task='tblastn')
                else:
                    target_short = '\nError while blasting with {}\n', format(index[n])
                short += target_short
                message += short + '  \n *** \n '

            st.success('Blast completed.')
            st.balloons()
            st.markdown(message)
            if path.exists(out_address + "{}_SBL.xml".format(job_name)):
                remove(out_address + "{}_SBL.xml".format(job_name))
            if path.exists(out_address + "{}.fasta".format(job_name)):
                remove(out_address + "{}.fasta".format(job_name))
