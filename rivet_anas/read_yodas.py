"""
usage: python3 read_yoda.py dir/filename [opt.][path [opt.]histogramnames]
"""

import yoda
import sys
from os.path import exists
import json

def read_histogram(file, histname, path='/ATLAS_2020_I1803608:CUTS=YES:TYPE=EW_ONLY/'):
    try: 
        hist = file[path+histname]
    except KeyError:
        print(f"No such histogram with path {path+histname}")
        return -1

    outdict = {}
    outdict["bins"] = list(hist.xEdges())
    outdict["yvals"] = list(hist.yVals())
    
    outfilename = '/data/horse/ws/maml087d-workspace/data_cuts/' + histname + '.json'
    i = 1

    # while exists(outfilename):
    #    outfilename =  '/data/horse/ws/maml087d-workspace/data_all' + str(i) + '/' + histname + ".json"
    #    i+=1

    with open(outfilename, 'w') as f:
        json.dump(outdict, f)

def main():
    if len(sys.argv) > 1:
        filename = str(sys.argv[1])
    else: 
        raise Exception("provide filename as commandline arg")

    if exists(filename):
        file = yoda.read(filename)
    else: raise Exception("file does not exits ... check path and filename")

    # check for specific histogram names in command line
    if len(sys.argv) == 3:
        path = str(sys.argv[2])
        namelist = []
        for i in range(9):
            namelist.append('temps_cos' +str(i))
            namelist.append('temps_phi' +str(i))

        for name in ['costheta', 'phi', 'xsec_bin'] + namelist:
            read_histogram(file, name, path=path)
        
    elif len(sys.argv) > 3:
        path = sys.argv[2]
        for name in sys.argv[3:]:
            print(name)
            read_histogram(file, name, path=path)
    else:
        # namelist = []
        # for i in range(9):
            # namelist.append('temps_cos' +str(i))
            # namelist.append('temps_phi' +str(i))

        for name in ['costheta', 'phi', 'xsec_bin'] #+ namelist:
            read_histogram(file, name)



if __name__ == "__main__":
    main()

