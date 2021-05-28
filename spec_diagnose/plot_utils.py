"""
docstring of plot_utils.
"""


import re

def AnnotateSegments(ax, RunDict, y=0., TerminationReason=False,
                     tref=0.):
    """
Annote a plot with information where segments begin.

  ax      - axes object to add annotations into
  RunDict - dictionary of data of run.  Created with ImportRun.
            must contain 'segs', 'tstart', 'termination', i.e.
            the data returned by FindLatestSegments

  y       - lower-bound y-value at which to print the (vertical) text objects
  TerminationReason
          - if ==True, print the termination reason of all segments which
            did not end with 'WallClock'
  tref    - use t-tref on x-axis (to avoid excessively long numbers)
"""

    segs=RunDict['segs']
    tstart=RunDict['tstart']
    termination=RunDict['termination']

    for seg,t,r in zip(segs,tstart,termination):
        ax.axvline(t-tref,color='grey',lw=0.5)
        # extract 'LevN_xx' from seg
        m=re.match('.*/(Lev._..)/Run',seg)
        if m:
            lev=m.group(1)
        else:
            print("WARNING.  segs has unexpected format")
            lev=""
        if TerminationReason and r not in ['WallClock',]:
            lev=lev+" "+r
        ax.text(t-tref, y,lev,rotation='vertical',verticalalignment='bottom', clip_on=True)


def PlotTruncationErrorSubdomain(ax, AdjustGrid, SD, tref=0., PileUpModes=False):
    """
Plot quantities relevant to assess truncation error for one subdomain.

ax -- axis object into which to plot data
AdjustGrid -- AdjustGrid dictionary, to be indexed by 'SD'
SD         -- name of the spherical shell to be plotted
tref       -- use t-tref as xaxis
"""

    a=AdjustGrid[SD] # shortcut

    # ==== get colors ====
    colors=[next(ax._get_lines.prop_cycler)['color'],
            next(ax._get_lines.prop_cycler)['color'],
            next(ax._get_lines.prop_cycler)['color']]

    # ==== construct labels ====
    labels=['0','1','2']
    if SD.startswith('Sphere'):
        labels=['r','theta','phi']
    if SD.startswith('FilledCylinder') or SD.startswith('Cylinder'):
        labels=['z','rho','phi']
        # see B2.hpp and DomainDetails.hpp, line ~ 67

    # plot extents
    t=a['Extents']['Extent[0]'][:,0]
    idx=0
    for key in a['Extents'].keys():
        if key=='time': continue
        N=a['Extents'][key][:,1]
        # 'pre':  at time-steps of adjustment, the value reported in
        #         AdjustGridDiagnostics.h5 is the *old* one.
        ax.step(t-tref, (N-20.)/10, where='pre', color=colors[idx],
                linewidth=2.5, label="(N{}-20)/10".format(labels[idx]))
        idx=idx+1


    # ==== plot basis-functions ====

    # Step 1: construct list of basisfunctions in correct order
    bfs=[bf for bf in a.keys() if bf.startswith('Bf')]
    bfs.sort()  # Harald thinks alphabetical sort is correct ...
    # ... except for B2Radial ...
    for k in range(len(bfs)):
        if 'B2Radial' in bfs[k]:
            if k!=2 or bfs[k]!='Bf1B2Radial' or bfs[1]!='Bf1B2':
                print("ERROR names/orders of basisfunctions unexpected. "
                      "Please check and amend this function."
                      "SD={}, bfs={}, k={}".format(SD,bfs,k))
                return
            # switch order (hardcode strings ok, as we only know how
            # to do this for these precise strings, cf. test just
            # above)
            bfs[1]='Bf1B2Radial'
            bfs[2]='Bf1B2'

    # Step 2: plot
    for idx,bf in enumerate(bfs):
        tmp=a[bf]['TruncationErrorExcess']
        ax.plot(tmp[:,0]-tref,tmp[:,1],'--',color=colors[idx],label='TruncErrExcess-{}'.format(labels[idx]),linewidth=1.5)
        tmp_idx=tmp[:,1]>0
        if sum(tmp_idx)>0:
            ax.plot(tmp[tmp_idx,0]-tref,tmp[tmp_idx,1],'o',color=colors[idx])
        if PileUpModes:
            tmp=a[bf]['MinNumberOfPiledUpModes']
            ax.plot(tmp[:,0]-tref,tmp[:,1],':', linewidth=1.5, color=colors[idx], label='# PileUpModes-{}'.format(labels[idx]))
    if tref==0:
        ax.set_xlabel('t/M')
    else:
        ax.set_xlabel(f'(t-{tref})/M')
    ax.legend();
    ax.set_title(SD)
    return


def PlotSubdomainConstraints(ax, GhCe, N=5, Ngrey=0):
    """
Make a plot of constraints.
  ax -- axes to plot into
  GhCe -- a dictionary containing constraint info, as obtained with segment_utils.LoadDat_from_segments
  title -- title of plot
  N -- plot the N subdomains with largest GhCe
  Ngrey -- plot the next 'Ngrey' subdomains in grey
"""
    maxD={}
    for legend,data in GhCe.items():
        if legend=='time': continue
        maxD[max(data[:,1])]=legend
    biggest=sorted(maxD.keys(),reverse=True)
    if N>len(biggest): N=len(biggest)
    if N+Ngrey>len(biggest): Ngrey=len(biggest)-N
    for b in biggest[:N]:
        legend=maxD[b]
        data=GhCe[legend]
        ax.plot(data[:,0],data[:,1],label=legend)

    for b in biggest[N:N+Ngrey]:
        legend=maxD[b]
        data=GhCe[legend]
        ax.plot(data[:,0],data[:,1],color='grey', lw=0.5)

    ax.set_xlabel('t/M')
    ax.legend()
    ax.set_yscale('log')



def PlotAH(ax, AH_dat, NormalizeRadii=True, title=None):
    """PlotAH(ax, AH)
    plot useful information about an apparent horizon.
    ax - axis object into which to plot the data
    AH - a dictionary with columns from Ah?.dat, as read by LoadDat_from_segments
 """

    # Plot rmin and rmax, possibly normalized
    if NormalizeRadii:
        norm=AH_dat['sqrt(Area/16pi)'][:,1]
        label_postfix='/Mirr'
    else:
        norm=1.
        label_postfix=''

    for q in 'min', 'max':
        # get a color for both curves
        color=next(ax._get_lines.prop_cycler)['color']
        tmp=q+'(r)'
        d=AH_dat[tmp]
        ax.plot(d[:,0],d[:,1]/norm,color=color, label=tmp+label_postfix)
        tmp=q+'(|r^i-c^i|)'
        d=AH_dat[tmp]
        ax.plot(d[:,0],d[:,1]/norm, '--', color=color, label=tmp+label_postfix)

    # plot remaining quantities
    d=AH_dat['sqrt(Area/16pi)']
    ax.plot(d[:,0],d[:,1],label='Mirr')

    d=AH_dat['L_surface']
    ax.plot(d[:,0],d[:,1]/10, label='L_surface/10')
    d=AH_dat['NumIterations']
    ax.plot(d[:,0],d[:,1]/10, lw=0.5, color='grey', label='Niterations/10')
    d=AH_dat['convg reason']
    ax.plot(d[:,0],d[:,1]/10, 'k--', lw=0.5, label='convg reason/10')

    ax.set_xlabel('t/M')
    ax.legend();

    if title is not None:
        ax.set_title(title,fontsize='x-large')


def PlotGravitationalWave(ax, waveform, l, m, label=None, title=None, RIndex=-1):
    """
Make a plot of gravitational waves.
  ax       -- axes to plot into
  waveform -- waveform data to plot (imported rh* or rPsi4* file with radii-entries)
  l,m      -- which mode to plot
  label    -- label for legend and axis, e.g. Psi4, h, M*Psi4, h/M
              Keep current label if option is not provided.
  RIndex -- index of wave extraction radius, default: -1 (outermost radius)

  Example:
    PlotGravitionalWave(ax, run['Psi4'], 2, 2, label='Psi4', RIndex=-1)
"""
    Ylm_str='Y_l'+format(l)+'_m'+format(m);
    Radius=list(waveform)[RIndex]
    Ylm=waveform[Radius][Ylm_str]

    # plot Re and Im over time
    keylist=list(Ylm)
    ReData=Ylm[keylist[1]]
    ImData=Ylm[keylist[2]]

    if label is not None:
        ax.plot(ReData[:,0],ReData[:,1], label='Re '+label+' Y'+format(l)+format(m)+'('+Radius+')')
        ax.plot(ImData[:,0],ImData[:,1], label='Im '+label+' Y'+format(l)+format(m)+'('+Radius+')')
        ax.legend()
        ax.set_ylabel(label)

    ax.set_xlabel('t/M')
    if title is not None:
        ax.set_title(title,fontsize='x-large')
