import os
import sys
import glob
import re
import h5py
import numpy as np


def FindLatestSegments(EvDir, Lev, interval=1e10):
    """
Search through a usual SpEC segment structure, to find all
segments within a time 'interval' of the last
existing segment.

Procedure:
1. Assemble list of existing directories
  [ ${EvDir}/Lev${Lev}_*] + [${EvDir}/Lev${Lev}_Ringdown/Lev${Lev}_* ]

2. Go backward through the directories
   2a. Retrieve start-time for each segment from DIR/Run/RestartTimes.txt
   2c. Terminate if start-time of current segment is
          earlier than ${DurationCompletedSegments} of start-time of last segment

3. This routine does NOT yet look for joined segments, JLev${Lev}.  Should be simple to add
JLev{$Lev} at start of list assembled in step 1., as long as CombineSegments also has combined Restart.txt


RETURNS
  segments,     -- list of *all* segments considered     [tuple of strings]
  tstart,       -- start time of each segment considered [tuple of floats]
  term_reason   -- TerminationReasons for each segment   [tuple of strings]
"""

    if not os.path.isdir(EvDir):
        raise IOError("Directory {} does not exist".format(EvDir))

    # find segments under consideration
    tmp=os.path.join(EvDir,"Lev{}_Ringdown/Lev{}_*".format(Lev,Lev),'Run')
    all_segments=sorted(glob.glob(tmp),reverse=True)

    tmp=os.path.join(EvDir,"Lev{}_*".format(Lev),'Run')
    all_segments.extend(sorted(glob.glob(tmp),reverse=True))

    segments=[]
    tstart=[]
    term_reason=[]
    for seg in all_segments:
        tmp=os.path.join(seg,'RestartTimes.txt')
        if not os.path.exists(tmp):
            raise IOError("{} not found -- don't yet know how to handle this".format(tmp))
        RestartTimes=np.loadtxt(tmp,
                                ndmin=1   # so it always returns an array
                               )

        tmp=os.path.join(seg,'TerminationReason.txt')
        if os.path.exists(tmp):
            with open(tmp, 'r') as myfile:
                TerminationReason = myfile.readlines()[0]
                prefix='Termination condition '
                if TerminationReason.startswith(prefix):
                    TerminationReason=TerminationReason[len(prefix):-1]
        else:
            TerminationReason='ongoing'

        segments.append(seg)
        tstart.append(RestartTimes[0])
        term_reason.append(TerminationReason)
        if tstart[0]-tstart[-1] > interval:
            break

    segments.reverse()
    tstart.reverse()
    term_reason.reverse()
    return segments, tstart,term_reason



def LoadH5_from_segments(segments, filename, dataset_matches='',group_matches='',
                         verbose=False):
    """
Given a list of segments (incl. '/Run' directories),
check each one for a file 'filename', load that h5 file, concatenate data, and
provide it as recursive dictionary

OPTIONS:
   dataset_matches=[regex] -- only load data-sets matching the regex (e.g. 'Y_l2_m2')
                               Note:  root-datasets are always loaded
   group_matches=[regex]   -- only traverse groups matching this regex (e.g. 'R0200')
"""

    # helper function for iterative loading
    def LoadFromOpenH5(F, D,dataset_matches='',group_matches=''):
        """Given the open H5-file handle F
        (either representing teh file, or a group inside the file),
        load all .dat files and store them as elements in the
        directory D.  Recursively descend into each .dir,
        and add those as sub-dictionaries into D
        dataset_matches=[regex] -- only load data-sets matching the regex"""

        for k in F.keys():
            if k.endswith('.dat') and \
            (re.match(dataset_matches,F[k].name)
             or F.parent==F # always keep top-level .dat fields
            ):
                # remove extension for convenience
                field=k[:-4]
                tmp=F[k][()]
                if field in D:
                    D[field]=np.concatenate((D[field], tmp))
                else:
                    D[field]=tmp
            if k.endswith('.dir') and re.match(group_matches,F[k].name):
                field=k[:-4]
                if field not in D:
                    D[field]={}
                LoadFromOpenH5(F[k],D[field],
                               dataset_matches=dataset_matches,
                               group_matches=group_matches)
        return

    D={}
    n_files=0
    for seg in segments:
        f=os.path.join(seg,filename)
        if os.path.exists(f):
            F=h5py.File(f,'r')
            LoadFromOpenH5(F,D,
                           dataset_matches=dataset_matches,
                           group_matches=group_matches)
            F.close()
            n_files=n_files+1
    if verbose:
        print("{} -- file found in {} of {} segments".format(filename,
                                                             n_files,len(segments)))
    return D



def LoadDat_from_segments_simple(segments, filename):
    """
Given a list of segments (incl. '/Run' directories),
check each one for a file 'filename', load that dat file, concatenate data, and
provide as np.array

Note: This function requires that the .dat file in all segments has
same number of columns
    """
    out=None
    n_files=0
    for seg in segments:
        f=os.path.join(seg,filename)
        if os.path.exists(f):
            tmp=np.loadtxt(f)
            if len(tmp.shape)==2: # ?? not sure why
                if out is None:
                    out=tmp
                else:
                    out=np.concatenate((out,tmp))
            n_files+=1
    print("{} -- file found in {} of {} segments".format(filename,
                                                         n_files,len(segments)))
    return out


def LoadDat_with_legend(F):
    """
Load a .dat file with standard SpEC legend strings.
Results will be placed into a dictionary indexed by the
legend string.

F -- filename

RETURNS
  D -- dictionary
"""
    if not os.path.isfile(F):
        raise IOError("File {} not found".format(F))

    # parse legend portion of file, construct legend as dictionary 'keys'
    keys={}   # dictonary of keys:  int -> legend string
    p=re.compile("^# *\[([0-9]+)\] * = *(.+)\n")
    for line in open(F):
        m=p.match(line)
        if m is not None:
            #print("{}  -{}-".format(m.group(1), m.group(2)))
            keys[ int( m.group(1) )]  = m.group(2).strip()
    #print(keys)
    tmp=np.loadtxt(F,
                   ndmin=2 # so it always returns a 2-d array
                  )
    #print(tmp.shape)
    D={}
    for idx,legend in keys.items():
        D[legend]=tmp[:,[0,idx-1]]  # -1, since SpEC legends are 1-based
    return D


def LoadDat_from_segments(segments, filename, verbose=False):
    """
Given a list of segments (incl. '/Run' directories),
check each one for a file 'filename', load that dat file, concatenate data, and
provide as np.array

The data is returned as a dictionary with 2-column data-sets, one each for each column,
named by the legend's inside the .dat files.  This allows for .dat files with a changing
number of columns, like GhCe_Linf.dat

RETURNS
  D -- dictionary
"""
    D={}

    n_files=0
    for seg in segments:
        f=os.path.join(seg,filename)
        if os.path.exists(f):
            tmp=LoadDat_with_legend(f)
            for legend, data in tmp.items():
                if legend in D:
                    D[legend]=np.concatenate((D[legend], data))
                else:
                    D[legend]=data
            n_files+=1
    if verbose:
        print("{} -- file found in {} of {} segments".format(filename,
                                                             n_files,len(segments)))
    return D