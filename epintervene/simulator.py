import numpy as np
import sandbox


def run():
    print('Running EpIntervene')
    print('Saving sample output')
    sample_data = [1, 2, 3]
    np.savetxt('./../data/sample_output.txt', np.array(sample_data))
    sandbox.run()
