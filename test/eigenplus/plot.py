import os
join = os.path.join
exists = os.path.exists
expanduser = os.path.expanduser
import sys
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pylab as plt
import numpy as np
tr = np.transpose
import pandas as pd

import mangan4 as mn

def run():

    dname = "_fig_rabi"
    if(not exists(dname)):
        os.makedirs(dname)
    
    man = mn.Man(root="_con_rabi")
    ts = man["t",:]
    calc  = man["cr",:] + 1j*man["ci",:]
    exact = man["exact_cr",:] + 1j*man["exact_ci",:]

    plt.plot(ts, abs(calc[:,0])**2,  "-",   label="calc0")
    plt.plot(ts, abs(calc[:,1])**2,  "-",   label="calc1")
    plt.plot(ts, abs(exact[:,0])**2, "k--", label="exact0")
    plt.plot(ts, abs(exact[:,1])**2, "k-.", label="exact1")
    plt.legend()
    plt.savefig(join(dname, "prob.pdf"))
    plt.close()

if __name__=="__main__":
    run()
    
    
