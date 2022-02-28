from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import streamlit as st
from numpy import arange
from Basic_Functions.primer_picker_funtion import *

def convert_rec(x):
    if isinstance(x, list):
        return list(map(convert_rec, x))
    else:
        return int(x)

def draw(exon, promoter, primer_pair):
    features = []
    count = 0

    for a in exon:
        count += 1
        features.append(GraphicFeature(start=a[0], end=a[1], strand=+1, color="#ffd700",
                                       label="Exon"))
    count = 0
    for a in promoter:
        count += 1
        features.append(GraphicFeature(start=a[0], end=a[1], strand=+1, color="#ffcccc",
                                       label="Promoter"))

    for a in primer_pair:
        tmp = primer_pair[a]
        features.append(GraphicFeature(start=tmp[0][0], end=tmp[0][1], strand=+1, color="#ccccff",
                                       label=a + '_F'))
        features.append(GraphicFeature(start=tmp[1][0], end=tmp[1][1], strand=-1, color="#ccccff",
                                       label=a + '_R'))
    return features

@st.cache
def read_exon(in_put):
    ls = in_put.replace(' ', '').split(',')
    tmp = []
    for a in ls:
        tmp.append(a.split(':', 1))
    return convert_rec(tmp)

@st.cache
def read_promoter(in_put):
    ls = in_put.replace(' ', '').split(',')
    tmp = []
    for a in ls:
        tmp.append(a.split(':', 1))
    return convert_rec(tmp)

@st.cache
def read_primer(in_put):
    pass


def dna_viewer():
    with st.form('main'):
        sample_seq = '''>sample_seq_1
            ggggtcttcccccaagaattggtagatgtggttttggtgccggccgacatgtctatagcttgtgagctctcgaatcccaatttaaagaagacaaaagcagaaaaacggatacca
            accaagattctgcttatccgggttttatgtggcttagttgttctctggctctgcttaagcttaagcttaggcttcttatgtatatgcaagaagaaagaggctgctgctgctgctga
            tgactcttcttcttctgctaaagggatgttgttcaggaatcagagcagaagtgagattgatgctatgctttctctcttctttgattcaaatcaggtttcccactacttctttcttc
            tttttattattgtttgacaatttaatttacttaaagcttttgaataattggagaacaggtaacatcttttgaatgtcgcaaggagaatggtggtattacatgttccttgtcaacac
            gttccgagaaaggagacgaggaggaggaggaggctaagagacatgttgttgcagaggttatgtcatcatctgagaatgaagaagaaggaggtgtgctgcatcaggttgtgttgttt
            tatgtaatgaacaaatggcattggtggttggtcctttgtgtactactagtgggcggcggccgtgtgatatttgtaagaaaagaagtttcttctttagtacaagacaagcagcagca
            gcagcagcaatgcaaaacagctgggaagtggaggaagaacatgcttctactcggcatcatcgcgggagtttctttgtctgttttatggttttgggataccaacgagaagatcttgt
            tccaaaggaaagagacgttaaccaacatgtgtgaggaacgagctcgggtgttgcaggaccagttcaatgttagcatgaaccatgtccacgccttgtccatcctcgtctctaccttt
            caccatggcaaaaccccttctgccattgatcaggtaaactaagagtccgaatggtaacatgtggtgcggttctaatagtaacaaaatcgtatagatacatgtgtaaatctagagat
            tttttaatgttaattgcggtgcgattctaatagtaaaaaaaagtatagatacatgtgtatatctagagatttttgttattgttaactacgatggtgttttaatagtaacaaaaacg
            tatagatacatatgtatatttagagatttttgttaatgttaattgcggcgcggttctaatagtaaccaaaaagtatatatacatgtgtatatctagagatttttttaactgcggtg
            ctgttccaatagtaacaaaaacgtatagatacatgtgtaaatctagagatttttgttaatgttaattgtggtgcagttctaatagtaacaaaaacgtatagatacatgtgtatatc
            taaagatttttgttaatgttaattgtggtgcggttccaatagtaacaaaaacgtatagatacatgtgtaaatctagagatttttgttaatgttaattgtggtgcggttctatagta
            acaaaaacgtatagatacatgtgtatatctagagatttttgttaatgctaattgtgctgcggtgcggttagaaataaaattccattcttaagatttttaactgttgattttttttt
            tttcag'''
        fasta = st.text_area('Input fasta sequence.', sample_seq, height=250)

        raw_seq, index = readFASTA(fasta)
        target_seq = seq_clean(raw_seq[0])
        # input exon,promoter,primer_pair
        exon_in = st.text_input('Please enter exon area:', value='300 : 600 , 900 : 1300', key='exon')
        promoter_in = st.text_input('Please enter promoter area:', value='10 : 250', key='promoter')
        primer_pair_in = st.text_input('Please enter promoter area:',
                                       value='200:220, 380:400,900:920,1080:1100',
                                       key='primer')
        primer_pair = {'pair1': [[200, 220], [380, 400]], 'pair2': [[900, 920], [1080, 1100]]}

        if st.form_submit_button('Draw gene'):
            exon = read_exon(exon_in)
            promoter = read_promoter(promoter_in)
            # primer_pair=read_primer(primer_pair_in)
            fig, (ax1, ax2) = plt.subplots(
                2, 1, figsize=(20, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]}
            )
            features = draw(exon, promoter, primer_pair)
            record = GraphicRecord(sequence_length=1600, features=features)
            record.plot(figure_width=20, ax=ax1)
            # PLOT THE LOCAL GC CONTENT (we use 50bp windows)
            record = Seq(target_seq)
            gc = lambda s: 100.0 * len([c for c in s if c in "GC"]) / 50
            xx = arange(len(record) - 50)
            yy = [gc(record[x: x + 50]) for x in xx]
            ax2.fill_between(xx + 25, yy, alpha=0.3)
            ax2.set_ylim(bottom=0)
            ax2.set_ylabel("GC(%)")
            st.pyplot(fig)
