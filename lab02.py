from IdFIR import SystemIdentification
import h5py
import numpy as np
import scipy.io

def main():
    data = sysID.load_file()
    nb = 10

    sysID.IdFIR(data, nb)


if __name__ == '__main__':
    sysID = SystemIdentification()
    main()