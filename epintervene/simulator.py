import numpy as np

#TODO: can remove this, and remove it from being included in the package

def run(fname):
    print('Running EpIntervene')
    print('Saving sample output')
    sample_data = [1, 2, 3]
    np.savetxt(fname, np.array(sample_data))
