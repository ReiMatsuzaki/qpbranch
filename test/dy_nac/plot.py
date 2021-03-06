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

    for name in ["tully1", "tully2"]:
        dname = "_fig_{0}".format(name)
        if(not exists(dname)):
            os.makedirs(dname)
    
        man = mn.Man(root="_con_{0}".format(name));
        t = man['t',:]
        q0 = man['q0',:]
        p0 = man['p0',:]
        c = man['c_re',:] + 1j*man['c_im',:]
        prob = man['prob',:]
        numI = man['_numI',0]
        
        for I in range(numI):
            plt.plot(t, prob[:,I], label=str(I))
        plt.plot(t, prob[:,0]+prob[:,1], "k--", label="total")
        plt.legend()
        plt.xlabel(r"$t$/fs", fontsize=15)
        plt.legend()
        plt.savefig(join(dname, "prob.pdf"))
        plt.close()
        
        plt.plot(t, q0)
        plt.xlabel(r"$t$/fs", fontsize=15)
        plt.savefig(join(dname, "q0.pdf"))
        plt.close()
        
        plt.plot(t, p0)
        plt.xlabel(r"$t$/fs", fontsize=15)
        plt.savefig(join(dname, "p0.pdf"))
        plt.close()

    man = mn.Man(root="_con_tully1".format(name));
    xs = man['xs',0]
    for v in ["v11", "v12", "v21", "v22"]:
        ys = man[v,0]
        plt.plot(xs, ys, label=v)
    plt.legend()
    plt.savefig("_fig_tully1/v.pdf")
    plt.close()

def run_delta():
    man = mn.Man(root="_con_delta_tully1")
    dir_fig = "_fig_delta_tully1"
    if(not exists(dir_fig)):
        os.makedirs(dir_fig)
        
    t = man['t',:]
    
    prob = man['prob',:]
    plt.plot(t, prob[:,0])
    plt.plot(t, prob[:,1])
    plt.savefig(join(dir_fig, "prob.pdf"))
    plt.close()

    q0 = man["q0",:]
    plt.plot(t, q0)
    plt.savefig(join(dir_fig, "q0.pdf"))
    plt.close()

    p0 = man["p0",:]
    plt.plot(t, p0)
    plt.savefig(join(dir_fig, "p0.pdf"))
    plt.close()    
    
def run_2state():
    man = mn.Man(root="_con_2state")
    ts = man["t",:]
    calc  = man["c_re",:] + 1j*man["c_im",:]
    exact = man["exact_cr",:] + 1j*man["exact_ci",:]

    dname = "_fig_2state"
    if(not exists(dname)):
        os.makedirs(dname)
        
    pl0, = plt.plot(ts, calc[:,0].real,  "-",   label="calc0")
    pl1, = plt.plot(ts, calc[:,1].real,  "-",   label="calc1")
    plt.plot(       ts, exact[:,0].real, "k--", label="exact0")
    plt.plot(       ts, exact[:,1].real, "k-.", label="exact1")
    plt.legend()
    plt.savefig(join(dname, "c_re.pdf"))
    plt.close()

    pl0, = plt.plot(ts, calc[:,0].imag,  "-",   label="calc0")
    pl1, = plt.plot(ts, calc[:,1].imag,  "-",   label="calc1")
    plt.plot(       ts, exact[:,0].imag, "k--", label="exact0")
    plt.plot(       ts, exact[:,1].imag, "k-.", label="exact1")
    plt.legend()
    plt.savefig(join(dname, "c_im.pdf"))
    plt.close()

    pl0, = plt.plot(ts, abs(calc[:,0])**2,  "-",   label="calc0")
    pl1, = plt.plot(ts, abs(calc[:,1])**2,  "-",   label="calc1")
    plt.plot(       ts, abs(exact[:,0])**2, "k--", label="exact0")
    plt.plot(       ts, abs(exact[:,1])**2, "k-.", label="exact1")
    plt.legend()
    plt.savefig(join(dname, "prob.pdf"))
    plt.close()

if __name__=="__main__":
    run()
    run_2state()
    run_delta()
