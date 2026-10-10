"""Portfolio explorer for saved gene-panel results; no private data required."""
import json
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
st.set_page_config(page_title='Less measurement, more focus', page_icon='ðŸ§¬', layout='wide')
report = json.loads((ROOT/'data/scanb_holdout_report.json').read_text())
st.title('How many genes are enough?')
st.markdown('A personal data science project exploring **selective gene measurement** for breast cancer subtype classification.')
st.caption('Explore saved research results. This demo does not diagnose patients or run a clinical PAM50 assay.')
overview, workflow, experiments, genes = st.tabs(['The result', 'How escalation works', 'The experiments', 'The gene panel'])
with overview:
    st.subheader('Start small. Measure more for uncertain cases.')
    c1,c2,c3 = st.columns(3)
    c1.metric('Staged balanced accuracy', f"{report['models']['staged']['balanced_accuracy']:.1%}")
    c2.metric('Average genes per patient', f"{report['average_genes']:.1f}", '47.2% fewer than 50')
    c3.metric('Patients needing the full panel', f"{report['escalation_fraction']:.1%}")
    st.write('On 916 reserved SCAN-B patients, the staged strategy used 20 genes first and expanded to 50 when its leading predictions were close. Measuring all 50 for everyone scored 93.2%; the staged approach scored 91.9%.')
    st.image(str(ROOT/'figures/scanb_staged_holdout.png'))
    with st.expander('What does balanced accuracy mean?'):
        st.write('Calculate the fraction correctly classified in each subtype, then average the four fractions. This gives smaller subtype groups the same importance as larger groups.')
    with st.expander('Limits of the result'):
        st.write('The observed 1.28-point loss met the planned 2-point target. Its 95% uncertainty interval was a loss of 0.37â€“2.33 points, so retention within 2 points is not established. Models were trained within SCAN-B: this is a separate-cohort replication, not direct transfer of TCGA models. Gene counts simulate measurement burden, not laboratory cost. Labels are PAM50-derived.')
with workflow:
    st.subheader('An illustrative routing example')
    st.write('The first model compares four subtype scores. If the gap between its two leading scores is below 0.45, the recorded strategy requests the remaining 30 genes.')
    gap=st.slider('Gap between the two leading scores',0.0,1.0,.20,.01)
    if gap < report['threshold']:
        st.info('Close scores â†’ add 30 genes â†’ use the 50-gene modelâ€™s prediction.')
    else:
        st.success('Wider score gap â†’ retain the 20-gene modelâ€™s prediction.')
    st.caption('This slider explains the frozen routing rule. It is not a patient prediction or a probability of being correct.')
    st.subheader('Why count only the additional genes?')
    fraction=st.slider('Illustrative fraction needing escalation',0,100,21)
    st.metric('Average distinct genes measured',f'{20+30*fraction/100:.1f}')
    st.caption('20 + 30 Ã— escalation fraction. Every first-stage gene is included in the larger panel.')
    st.write('On the actual holdout, extra measurements corrected 42 predictions, introduced 16 errors, and left 29 initially incorrect escalated predictions incorrect.')
with experiments:
    st.subheader('The path to the final experiment')
    st.dataframe(pd.DataFrame([
        {'Experiment':'TCGA automatic 20 â†’ 100','Average genes':43.5,'Staged accuracy':'86.3%','Full panel':'88.5%','Observed loss':'2.24 points','Evaluation':'Nested development'},
        {'Experiment':'TCGA PAM50 subset 20 â†’ 50','Average genes':22.2,'Staged accuracy':'90.8%','Full panel':'93.1%','Observed loss':'2.33 points','Evaluation':'Nested development'},
        {'Experiment':'SCAN-B PAM50 subset 20 â†’ 50','Average genes':26.4,'Staged accuracy':'91.9%','Full panel':'93.2%','Observed loss':'1.28 points','Evaluation':'Fresh 916-case holdout'},
    ]),hide_index=True,width='stretch')
    st.write('The first two designs missed the 2-point goal. Before SCAN-B modeling, the threshold-selection rule became more conservative: retain performance within 1 point during development, then assess the original 2-point target on reserved patients.')
    st.caption('These experiments differ in cohort and development rule; the comparison does not isolate which change caused the improvement.')
    st.subheader('Original gene-count comparison')
    st.image(str(ROOT/'figures/panel_size_curve.png'))
    st.caption('Training-validation scores. PAM50 markers represent models using its fixed gene list, not the original PAM50 classifier.')
with genes:
    st.subheader('The final 20-gene first stage')
    panel=pd.read_csv(ROOT/'data/scanb_final_small_panel.csv')
    st.dataframe(panel.rename(columns={'scanb_symbol':'Gene','feature_id':'Entrez feature'}),hide_index=True,width='stretch')
    st.download_button('Download the 20-gene list',panel.to_csv(index=False),'first_stage_genes.csv','text/csv')
    st.write('ANOVA ranked differences in expression between subtypes using development patients only. The full second stage uses all 50 PAM50 genes. The same first-stage list is used for every patient; escalation depends on model scores, not knowing the true subtype.')
st.divider()
st.caption('Python Â· scikit-learn Â· TCGA Â· SCAN-B Â· AI-assisted implementation. Source and reproduction steps are available in the repository.')
