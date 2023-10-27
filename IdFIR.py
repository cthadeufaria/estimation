import matlab.engine
import numpy as np
import h5py
import scipy.io
# from scipy.linalg import pinv
# from scipy.signal import lfilter


class SystemIdentification:
    def __init__(self) -> None:
        self.eng = matlab.engine.start_matlab()

    def IdFIR(self, data, nb):
        Y = data['y'][nb + 1:]
        Phi = np.zeros((len(data['u']) - nb, nb + 1))

        for i in range(len(data['u']) - nb):
            Phi[i, :] = data['u'][i + nb:i:-1]

        Theta = np.linalg.lstsq(Phi, Y, rcond=None)[0]
        # Alternatively, you can use the commented line below to calculate Theta using the normal equation
        # Theta = np.linalg.solve(Phi.T @ Phi, Phi.T @ Y)

        sys = {'num': [0] + Theta.tolist(), 'den': [1]}

        return Y, Phi, sys
    
    def load_file(self, path):
        try:
            data = self.eng.load(path, nargout=1)
        except:
            data = scipy.io.loadmat(path)
        else:
            file = h5py.File(path, 'r')
            data = file.get('data/variable1')
            data = np.array(data)

        return data