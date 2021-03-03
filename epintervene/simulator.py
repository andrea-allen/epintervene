import numpy as np


def run(fname):
    print('Running EpIntervene')
    print('Saving sample output')
    sample_data = [1, 2, 3]
    np.savetxt(fname, np.array(sample_data))
