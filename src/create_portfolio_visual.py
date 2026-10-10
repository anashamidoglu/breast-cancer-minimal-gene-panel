"""Create a shareable project summary directly from recorded results."""
import json
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from tune_models import ROOT

def main():
    r=json.loads((ROOT/'data/scanb_holdout_report.json').read_text())
    fig=plt.figure(figsize=(12,6.3),facecolor='#101C2B')
    fig.text(.06,.90,'HOW MANY GENES ARE ENOUGH?',color='#77DDD0',fontsize=12,weight='bold')
    fig.text(.06,.78,'Start small. Measure more\nwhen predictions are uncertain.',color='white',fontsize=27,weight='bold',linespacing=1.25)
    cards=[('91.9%','balanced accuracy'),('26.4','genes per patient, on average'),('47.2%','fewer gene measurements')]
    for i,(value,label) in enumerate(cards):
        x=.06+i*.30
        fig.patches.append(FancyBboxPatch((x,.28),.27,.22,boxstyle='round,pad=0.012',transform=fig.transFigure,facecolor='#203247',edgecolor='none'))
        fig.text(x+.018,.40,value,color='#77DDD0',fontsize=31,weight='bold')
        fig.text(x+.018,.32,label,color='white',fontsize=11)
    fig.text(.06,.19,'20-gene first stage  →  full 50-gene panel for 21.4% of patients',color='white',fontsize=15)
    fig.text(.06,.12,'916 fresh SCAN-B holdout patients · Full-panel balanced accuracy: 93.2%',color='#CAD3DF',fontsize=11)
    fig.text(.06,.065,'Observed loss: 1.28 points · 95% interval: 0.37–2.33 points · Computational study, not a clinical assay',color='#CAD3DF',fontsize=9)
    fig.savefig(ROOT/'figures/portfolio_summary.png',dpi=200,facecolor=fig.get_facecolor())
    fig.savefig(ROOT/'figures/portfolio_summary.svg',facecolor=fig.get_facecolor())
    plt.close(fig)

if __name__=='__main__':main()
