"""
read the angular coeffs from the output yoda
usage: python3 read_angcoeffs.py dir/yodafile outpath yodapath ptbinpath
"""

import yoda
import sys
from os.path import exists
import json


def read_ang(file, tprofile, path='/ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/', outpath='/data/horse/ws/maml087d-workspace/rivet_anas/data/coeffs.txt'):
    if path+tprofile in file.keys():
        tpf = file[path+tprofile]
        print(f"exporting {path+tprofile}", end=' ..... ')
    elif '/raw'+path+tprofile in file.keys():
        tpf = file['/raw'+path+tprofile]
        print(f"exporting {'/raw'+path+tprofile}", end=' ..... ')

    else:
        print(f"No such Tprofile with path {path+tprofile} or {'/raw'+path+tprofile}")
        return -1

    # if isinstance(tpf, yoda.core.Profile1D):
    coeff = tpf.yVals()
    with open(outpath, "a") as f:
        for c in coeff:
            f.write(str(c)+'\n')
    print("done")
    # else:
    #     print("not a Profile1D")
    #     return -2



def main():
    print(__doc__)

    if len(sys.argv) > 1:
        filename = str(sys.argv[1])
    else:
        raise Exception("provide filename as commandline arg")

    if len(sys.argv) > 2:
        outpath = str(sys.argv[2])
    else:
        raise Exception("provide path for outfiles")

    if exists(filename):
        file = yoda.read(filename)
    else: raise Exception("file does not exits ... check path and filename")

    #clearing file
    with open(outpath+"coeffs.txt", "w") as f:
        f.write("")
    with open(outpath+"coeffs_of.txt", "w") as f:
        f.write("")

    if len(sys.argv) == 5:
        path = sys.argv[3]
        ptbins = sys.argv[4]
        namelist = []

        for i in range(8):
            namelist.append(f"A_{i}")
        if exists(ptbins):
            with open(ptbins, 'r') as f:
                n = sum(1 for line in f)
        else: print(f"{ptbins} doesnt exist"); raise Exception()
        for i in range(n-1):
            namelist.append(f"datbin{i}_xsec")

        for name in namelist:
            read_ang(file, name, path=path, outpath=outpath+'coeffs.txt')

        for name in [f"overflow_ang{i}" for i in range(8)] + ["overflow_xsec"]:
            read_ang(file, name, path=path, outpath=outpath+'coeffs_of.txt')

    else:
        print("not implemented yet")






if __name__ == "__main__":
    main()
