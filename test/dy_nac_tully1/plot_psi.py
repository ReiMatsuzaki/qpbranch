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
    xs = man['xs',0]

    its = [0, 20, 40, 60, 80]
    for it in its:
        psi0 = man['psi_0_re',it] + 1j*man['psi_0_im',it]
        psi1 = man['psi_1_re',it] + 1j*man['psi_1_im',it]
        p, = plt.plot(xs, psi0.real, label="1, it={0}".format(it))
        plt.plot(     xs, psi1.real, linestyle="--", color=p.get_color(),
                 label="2, it={0}".format(it))
    plt.legend()
    plt.savefig("psi.pdf")

if __name__=="__main__":
    run()
    
