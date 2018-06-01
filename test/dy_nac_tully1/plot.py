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
    man = mn.Man();
    #    t = man['t',:]*mn.param.get()["u.FS"]
    t = man['t',:]
    q0 = man['q0',:]
    p0 = man['p0',:]
    c = man['c_re',:] + 1j * man['c_im',:]
    prob = man['prob',:]
    numI = man['_numI',0]

    for I in range(numI):
        plt.plot(t, prob[:,I], label=str(I))
    plt.plot(t, prob[:,0]+prob[:,1], "k--", label="total")
    plt.legend()
    plt.xlabel(r"$t$/fs", fontsize=15)
    plt.legend()
    plt.savefig("prob.pdf")
    plt.close()

    plt.plot(t, q0)
    plt.xlabel(r"$t$/fs", fontsize=15)
    plt.savefig("q0.pdf")
    plt.close()

    plt.plot(t, p0)
    plt.xlabel(r"$t$/fs", fontsize=15)
    plt.savefig("p0.pdf")
    plt.close()

if __name__=="__main__":
    run()
    
    
