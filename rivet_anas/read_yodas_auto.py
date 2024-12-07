"""
usage: python3 read_yoda_auto.py dir/filename outpath [opt.][path [opt.]histogramnames]
"""

import yoda
import sys
from os.path import exists
import json

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
    print(f"exporting histo: {outfilename}")
    # while exists(outfilename):
    #    outfilename =  '/data/horse/ws/maml087d-workspace/data_all' + str(i) + '/' + histname + ".json"
    #    i+=1

    with open(outfilename, 'w') as f:
        json.dump(outdict, f)

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
    if len(sys.argv) == 4:
        path = str(sys.argv[3])
        namelist = []
        
        for i in range(9):
            namelist.append('temps_cos' +str(i))
            namelist.append('temps_phi' +str(i))
            namelist.append('temps_2D' + str(i))
        print(outpath)
        for name in ['costheta', 'phi', 'xsec_bin'] + namelist:
            read_histogram(file, name, path=path, outpath=outpath)
        
    elif len(sys.argv) > 4:
        path = sys.argv[3]
        for name in sys.argv[4:]:
            print(name)
            read_histogram(file, name, path=path, outpath=outpath)
    else:
        namelist = []
        for i in range(9):
            namelist.append('temps_cos' +str(i))
            namelist.append('temps_phi' +str(i))

        for name in ['costheta', 'phi', 'xsec_bin'] + namelist:
            read_histogram(file, name, outpath=outpath)



if __name__ == "__main__":
    main()

