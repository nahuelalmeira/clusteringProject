import os
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', 'k', 'y',
          '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', 'k', 'y']
markers = ['s', 'o', 'v', '^', '*', '<', '>', '8', '.', '',
		   's', 'o', 'v', '^', '*', '<', '>', '8', '.', '']

rc_font_size = 30
rc_label_size = 28
rc_legend_size = 26
mpl.rcParams['figure.figsize'] = (12, 8)
mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['legend.fontsize'] = rc_legend_size
mpl.rcParams['savefig.transparent'] = True
mpl.rcParams['savefig.bbox'] = 'tight'
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['axes.labelsize'] = rc_font_size
mpl.rcParams['xtick.labelsize'] = rc_label_size
mpl.rcParams['ytick.labelsize'] = rc_label_size
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['xtick.minor.width'] = 2
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['ytick.minor.width'] = 2

mpl.rcParams['patch.linewidth'] = 1.5
mpl.rcParams['figure.titlesize'] = 24

#from matplotlib import rc
#font = {'family'     :  'sans-serif',
#        'sans-serif' : ['Helvetica'] }

#rc('font',**font)

def test1(savefig=False):

    x = np.linspace(0, 5, 100)
    y = np.sin(x)

    plt.figure()
    plt.xlabel("$x$")
    plt.ylabel("$\mathrm{sin}\,(x)$")
    plt.plot(x, y, label="$\mathrm{Trigonometric\;function}")
    plt.plot(x, y+1, label="Trigonometric function 2")
    plt.legend(loc="best")
    if savefig:
        plt.savefig("test1.pdf")
    else:
        plt.show()

if __name__ == '__main__':

    print(os.getcwd())
    test1()
    test1(savefig=True)
