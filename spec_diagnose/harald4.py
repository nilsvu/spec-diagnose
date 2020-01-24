from __future__ import print_function
import glob
import re
import os
import sys
import h5py
import numpy as np
from collections import namedtuple


def convert(obj,name='GenericDict',depth=1111):
    """convert a dictionary into a namedtuple
    Taken from https://gist.github.com/href/1319371
    """
    if isinstance(obj, dict):
        if depth>0:
            # convert entries
            for key, value in obj.items():
                obj[key] = convert(value,key,depth=depth-1)
        return namedtuple(name, obj.keys())(**obj)
    elif isinstance(obj, list):
        if depth>0:
            return [convert(item, name, depth=depth-1) for item in obj]
        else:
            return obj
    else:
        return obj

    
def FindFiles(path, Lev, filename):
    """identify all .h5 files to be loaded by this LoadH5 command.
Procedure:
1a. Look for Lev{Lev}_* and find the files there
1b. Look for Lev{Lev}_Ringdown/Lev{Lev}_* to find the files there
1c. Combine list of 2a and 2b

2. Look for the file in JLev{Lev}, which is supposed to contain 
   joined files
2a. If that file exists, replace those files from 1c which are **older** than 
    JLev{Lev}'s file by JLev{Lev}
2b. If that file does not exist, suggest to run CombineSegments

RETURNS  
  list_of_files, JLev_present
     list_of_files -- list of files to read
     JLev_present  -- bool to indicate whether JLev{} is used

"""
    
    if not os.path.isdir(path):
        raise IOError("Directory {} does not exist".format(path))

    # find files under consideration
    tmp=os.path.join(path,"Lev{}_*".format(Lev),'Run',filename)
    Lev_files=sorted(glob.glob(tmp))
    tmp=os.path.join(path,"Lev{}_Ringdown/Lev{}_*".format(Lev,Lev),'Run',filename)
    Lev_files.extend(sorted(glob.glob(tmp)))
    #print "files={}".format(files)
    if len(Lev_files)==0:
        raise IOError("No files found that match \'{}\'".format(filename))

    jfile=os.path.join(path,"JLev{}".format(Lev),filename)
    if os.path.isfile(jfile):
        #print(">> Using JLev{}/{} for speed-up".format(Lev,filename))
        age_J=os.path.getmtime(jfile)
        new_Lev_files=[f for f in Lev_files if os.path.getmtime(f) > age_J]
        out=[jfile]+new_Lev_files
        JLev_used=True
    else:
        out=Lev_files
        JLev_used=False
    if len(out)>20:
        print("| {} Lev{} is big.  Run".format(path,Lev))
        print("| $ CombineSegments.py -o JLev{Lev} -L{Lev} --gnuparallel -x 'Profiler|PowerDiag|FilterDiag|CceR.....h5|Phi'".format(Lev=Lev))
    return out,JLev_used


def LoadH5(path, Lev, filename, dataset_matches='',group_matches=''):
    """Search segments for files name 'filename'
    as given by path & Lev.  Concatenate data in 'filename' in all segments, 
    and provide it as a recursive dictionary.

    OPTIONS: 
       dataset_matches=[regex] -- only load data-sets matching the regex (e.g. 'Y_l2_m2')
                                Note:  root-datasets are always loaded
       group_matches=[regex]   -- only traverse groups matching this regex (e.g. 'R0200')
  
Convenience features:
      - Run/Lev{Lev}  inserted automatically
      - for 'Horizons.h5',  'ApparentHorizons/' inserted automatically
      - filename='rh' is expanded to 'GW2/rh_FiniteRadii_CodeUnits.h5'
      - filename='rPsi4' is expanded to 'GW2/rPsi4_FiniteRadii_CodeUnits.h5'
      - filename='AdjustGridExtents.h5'.  Always return as dictionary (b/c too big for namedtuple)
    """

    if filename in ['Horizons.h5', 'RedshiftQuantities.h5']:
        filename=os.path.join('ApparentHorizons', filename)
    if filename=='rh':
        filename='GW2/rh_FiniteRadii_CodeUnits.h5'
    if filename=='rPsi4':
        filename='GW2/rPsi4_FiniteRadii_CodeUnits.h5'
    if not os.path.isdir(path):
        raise IOError("Directory {} does not exist".format(path))

        
    # find files under consideration
    #tmp=os.path.join(path,"Lev{}_*".format(Lev),'Run',filename)
    #files=sorted(glob.glob(tmp))
    #tmp=os.path.join(path,"Lev{}_Ringdown/Lev{}_*".format(Lev,Lev),'Run',filename)
    #filesRD=sorted(glob.glob(tmp))
    #files.extend(filesRD)
    #print "files={}".format(files)
    #if len(files)==0:
    #    raise IOError("No files found that match \'{}\'".format(filename))

    files, JLev=FindFiles(path,Lev,filename)
    
    # helper function for iterative loading
    def LoadFromOpenH5(F, D,dataset_matches='',group_matches=''):
        """Given the open H5-file handle F
        (either representing / or a group inside the file),
        load all .dat files and store them as elements in the
        directory D.  Recursively descend into each .dir -group,
        and add those as sub-dictionaries into D

        dataset_matches=[regex] -- only load data-sets matching the regex"""

        #print "A-- F.keys()=",F.keys()
        for k in F.keys():
            #print k
            #print(F[k].name)
            if k.endswith('.dat') and \
            (re.match(dataset_matches,F[k].name)
             or F.parent==F # always keep top-level .dat fields
            ):
                # remove extension, since periods cannot be
                # used in tuple-fields.  Also, convenient!
                field=k[:-4]
                # also replace '-' by 'm', since namedtuple
                # doesn't do -.  This happens in Y_l2_m-2.dat
                #field=field.replace('-','m')
                #field=field.replace('.','x')
                tmp=F[k][()]
                if field in D:
                    idx=np.argmax(tmp[:,0]>D[field][-1,0])
                    D[field]=np.concatenate((D[field], tmp[idx:,:]))
                else:
                    D[field]=tmp
            if k.endswith('.dir') and re.match(group_matches,F[k].name):
                field=k[:-4]
                #field=field.replace('.','x')
                if field not in D:
                    D[field]={}
                LoadFromOpenH5(F[k],D[field],
                               dataset_matches=dataset_matches,
                               group_matches=group_matches)
        return


    # first fill a dictionary, which is extensible
    D={}
    Nfiles=len(files)
    width=max(20-Nfiles, 0)
    n=0
    print("."*Nfiles,end='')
    Jstring="J" if JLev else "="
    for f in files:
        n=n+1
        print('\r'+Jstring+"="*(n-1)+"."*(Nfiles-n)+" "*width+" {}".format(os.path.relpath(f,path)),
              end='')
  
        F=h5py.File(f,'r')
        #print F.keys()
        LoadFromOpenH5(F,D,
                       dataset_matches=dataset_matches,
                       group_matches=group_matches)
        F.close()
    print()
    return D


def LoadDat(path, Lev, filename, dataset_matches='',group_matches=''):
    """Search segments for files name 'filename'
    as given by path & Lev.  Concatenate data in 'filename' in all segments, and return it.
  
    """


    if not os.path.isdir(path):
        raise IOError("Directory {} does not exist".format(path))

    
    files, JLev=FindFiles(path,Lev,filename)
    
    # first fill a dictionary, which is extensible
    out=None
    Nfiles=len(files)
    width=max(20-Nfiles, 0)
    n=0
    print("."*Nfiles,end='')
    Jstring="J" if JLev else "="
    for f in files:
        n=n+1
        print('\r'+Jstring+"="*(n-1)+"."*(Nfiles-n)+" "*width+" {}".format(os.path.relpath(f,path)),end='')
        #print(f)
        tmp=np.loadtxt(f)
        #print(tmp.shape)
        #print(len(tmp.shape))
        #print(tmp)
        if len(tmp.shape)==2:
            if out is None:
                out=tmp
            else:
                #print
                #print("out.shape={}, tmp.shape{}".format(out.shape, tmp.shape))
                #print
                idx=np.argmax(tmp[:,0]>out[-1,0])
                out=np.concatenate((out,tmp[idx:,:]))
    print()
    return out

