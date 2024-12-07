"""
usage: python3 read_yoda.py dir/filename outpath yodapath ptbinpath
"""

import yoda
import sys
from os.path import exists, isdir
from os import mkdir
import json
# import ipdb

def jsondump(data, file):
    json.dump(data, file, indent=4)

def read_histogram(file, histname, outpath='/data/horse/ws/maml087d-workspace/data_cuts', path='/ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/'):
    if path+histname in file.keys():
        hist = file[path+histname]
    elif '/raw'+path+histname in file.keys():
        hist = file['/raw'+path+histname]
    else:
        print(f"No such histogram with path {path+histname} {'/raw'+path+histname}")
        return -1
        
#    try: 
#        hist = file[path+histname]
#    except KeyError:
#        print(f"No such histogram with path {path+histname}")
#        return -1

    outdict = {}
    if isinstance(hist, yoda.core.Histo1D):
        outdict["bins"] = list(hist.xEdges())
        outdict["yvals"] = list(hist.yVals())
    elif isinstance(hist, yoda.core.Histo2D):
        outdict["bins"] = [list(hist.xEdges()), list(hist.yEdges())]
        outdict["zvals"] = list(hist.zVals())
    else: print("not a 1DHisto or 2DHisto")
    outfilename = outpath + '/' +histname + '.json'
    i = 1
    print(f"exporting histo: {outfilename} ......", end='')
    # while exists(outfilename):
    #    outfilename =  '/data/horse/ws/maml087d-workspace/data_all' + str(i) + '/' + histname + ".json"
    #    i+=1

    with open(outfilename, 'w') as f:
        jsondump(outdict, f)

    print("done")


def read_binned_temps(file, ptbinpath='/data/horse/ws/maml087d-workspace/rivet_anas/data/ptbins.txt', path='/ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/', outpath=None):

    for j in range(9):
        for i in range(countlines(ptbinpath)-1):
            for k in ["_cos", "_phi", "_2D"]:
                read_histogram(file, "ptbin"+str(i)+k+str(j), outpath=outpath, path=path)
                read_histogram(file, "ptoverflow"+k+str(j), outpath=outpath, path=path)


def read_data(file,  ptbinpath='/data/horse/ws/maml087d-workspace/rivet_anas/data/ptbins.txt', path='/ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/', outpath=None):
    print(f"checking if histofolder {outpath} exists")
    if not exists(outpath):
        s = str(input("It doesnt exist ... do u want me to make the dir? [y|n]")).lower()[0]
        if s == 'y':
                makedirs(outpath)
        else:
            print("pls use an existing folder as outpath + histofolder")
            return
    else:
        if not isdir(outpath):
            print("{outpath} is not a directory")
            return

    for j in ["_cos", "_phi", "_xsec", "_2D"]:
        for i in range(countlines(ptbinpath)-1):
            read_histogram(file, "datbin"+str(i)+j, outpath=outpath, path=path)
        read_histogram(file, "overflow"+j, outpath=outpath, path=path)


def countlines(filepath='/data/horse/ws/maml087d-workspace/rivet_anas/data/ptbins.txt'):
    if exists(filepath):
            with open(filepath, 'r') as f:
                n = sum(1 for line in f)
    return n

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

    # check for specific histogram names in command line
    if len(sys.argv) == 5:
        path = str(sys.argv[3])
        ptbinpath = str(sys.argv[4])

        print(f"reading: {filename} with path {path} and ptbins={ptbinpath}")

        # ipdb.set_trace()
        for k in file.keys():
            if k.find("datbin") != -1:
                read_data(file, ptbinpath=ptbinpath, outpath=outpath, path=path)
                return 0
            elif k.find("ptbin") != -1:
                read_binned_temps(file, ptbinpath=ptbinpath, outpath=outpath, path=path)
                return 0
        else:
            print("no histo with name '*datbin*' or '*ptbin*' found -> terminating without reading yodas")
            return -1
    else:
        print("wrong number of args")
        print(__doc__)



if __name__ == "__main__":
    main()

