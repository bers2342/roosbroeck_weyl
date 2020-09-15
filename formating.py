import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter



def step1(fontsize=22,pdf=3) :
    mpl.rcdefaults()

    mpl.rcParams['mathtext.default']= 'regular'

    mpl.rcParams['font.size'] = fontsize # change the size of the font in every figure

    mpl.rcParams['font.family'] = 'serif' # font Arial in every figure

    mpl.rcParams['font.weight'] = 100 # font Arial in every figure

    mpl.rcParams['axes.labelsize'] = fontsize

    mpl.rcParams['xtick.labelsize'] = fontsize

    mpl.rcParams['ytick.labelsize'] = fontsize

    mpl.rcParams['axes.linewidth'] = 0.6 # thickness of the axes lines

    mpl.rcParams['pdf.fonttype'] = pdf  # Output Type 3 (Type3) or Type 42 (TrueType), TrueType allows

                                    # editing the text in illustrator

def set_ax(ax,xmin,xmax,ymin,ymax,xlabel,ylabel,title="") :
    
    ax.set_xlim(xmin,xmax)   # limit for xaxis
    ax.set_ylim(ymin,ymax) # leave the ymax auto, but fix ymin
    ax.set_xlabel(xlabel, labelpad = 5)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    #############################################
    ## Allow to shift the label ticks up or down with set_pad ##

    for tick in ax.xaxis.get_major_ticks():        
        tick.set_pad(7)

    for tick in ax.yaxis.get_major_ticks():
        tick.set_pad(8)

def set_ticks(ax,space_x,space_y,length_major=7,length_minor=3,xminor=True,yminor=True) :
    xtics = space_x # space between two ticks
    mxtics = xtics / 2  # space between two minor ticks
    
    ytics = space_y
    mytics = ytics / 2   # or "AutoMinorLocator(2)" if ytics is not fixed, just put 1 minor tick per interval

    majorFormattery = FormatStrFormatter('%g') # put the format of the number of ticks
    majorFormatterx = FormatStrFormatter('%g')
    
    ax.tick_params('both', length=length_major, width=1, which='major',direction='in',right='on',top='on')
    ax.tick_params('both', length=length_minor, width=1, which='minor',direction='in',right='on',top='on')

    ax.xaxis.set_major_locator(MultipleLocator(xtics))
    ax.xaxis.set_major_formatter(majorFormatterx)
    if xminor :
        ax.xaxis.set_minor_locator(MultipleLocator(mxtics))
        
    ax.yaxis.set_major_locator(MultipleLocator(ytics))
    ax.yaxis.set_major_formatter(majorFormattery)
    if yminor :
        ax.yaxis.set_minor_locator(MultipleLocator(mytics))