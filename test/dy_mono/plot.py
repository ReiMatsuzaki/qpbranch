import os
join = os.path.join
exists = os.path.exists
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pylab as plt

import mangan4 as mn

def run(name):
    man = mn.Man(root="_con_{0}".format(name))

    t = man["t",:]
    
    for lbl in ["q0", "p0", "gr0", "gi0"]:
        calc = man[lbl,:]
        exact = man["exact_"+lbl,:]
        plt.plot(t, calc,  "-",  label="calc")
        plt.plot(t, exact, "--", label="exact")
        plt.legend()
        dirname = "_fig_"+name
        if(not exists(dirname)):
            os.mkdir(dirname)
        plt.savefig(join(dirname, "{0}.pdf".format(lbl)))
        plt.close()

if __name__=="__main__":
    run("fp")
    run("harm")
    
