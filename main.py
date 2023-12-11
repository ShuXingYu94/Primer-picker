from PIL import Image
from primer_picker import *
from s_haplotype import *
import streamlit as st
from in_progress import *

APP_TITLE = "Primer Picker"

# Set the configs
st.set_page_config(
    page_title=APP_TITLE,
    # page_icon=Image.open(r'./primer_picker_cache/Icon/word_P.ico'),
    layout="wide",
    initial_sidebar_state="auto",
)


# icon = Image.open(r'./primer_picker_cache/Icon/black on white.png')
# report = get_system_report()

def general_main(icon):
    st.markdown(""" <style>
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style> """, unsafe_allow_html=True)
    st.title('Primer Picker and Specificity Check')
    st.sidebar.image(icon, use_column_width=True)
    st.sidebar.markdown('''[![ShuXingYu94 - Primer-picker](https://img.shields.io/static/v1?label=ShuXingYu94&message=Primer-picker&color=green&logo=github)](https://github.com/ShuXingYu94/Primer-picker "Go to GitHub repo")
        [![GitHub tag](https://img.shields.io/github/tag/ShuXingYu94/Primer-picker?include_prereleases=&sort=semver&color=green)](https://github.com/ShuXingYu94/Primer-picker/releases/)''')
    st.info('In order to work smoothly, please check the information on the left.')

# Main function
def Primer_Picker_Main():
    # Make general page information.
    general_main(icon)

    # Choose from different fucntionalities.
    state = st.sidebar.selectbox('Function', ['Primer Picker (Single)','S-haplotype Blast','DNA Viewer', 'Primer Picker (Multiple)', 'DNA to Amino Acid'
                                              ])
    if state == 'Primer Picker (Single)':
        primer_picker_single_main()
    elif state== 'S-haplotype Blast':
        s_haplotype_main()
    elif state == 'Primer Picker (Multiple)':
        st.sidebar.warning('Currently not available.')
        primer_picker_multiple_main()
    elif state == 'DNA to Amino Acid':
        st.sidebar.warning('Currently not available.')
        DNA_to_Amino()
    elif state == 'DNA Viewer':
        import dna_viewer
        dna_viewer.dna_viewer()


# Run the Primer_Picker
if __name__ == '__main__':
    try:
        Primer_Picker_Main()
    except (ValueError, IndexError) as val_ind_error:
        st.error(f"There is a problem with values/parameters or dataset due to {val_ind_error}.")
    except TypeError as e:
        # st.warning("TypeError exists in {}".format(e))
        pass
