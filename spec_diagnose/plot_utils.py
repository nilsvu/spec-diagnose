import re

def AnnotateSegments(ax, segs, tstart, termination, y=0.):
    """
Annote a plot with information where segments begin.

  ax - axes object to add annotations into
  segs
  tstart
  termination -- data returned by segment_utils.FindLatestSegments

  y - lower-bound y-value at which to print the (vertical) text objects
"""
    for seg,t in zip(segs,tstart):
        ax.axvline(t,color='grey',lw=0.5)
        # extract 'LevN_xx' from seg
        m=re.match('.*/(Lev._..)/Run',seg)
        if m:
            lev=m.group(1)
        else:
            print("WARNING.  segs has unexpected format")
            lev=""
        ax.text(t, y,lev,rotation='vertical',verticalalignment='bottom', clip_on=True)


def PlotTruncationError_Sphere(ax, AdjustGrid, SD):
    """
Plot quantities relevant to assess truncation error for one spherical shell.

ax -- axis object into which to plot data
AdjustGrid -- AdjustGrid dictionary, to be indexed by 'SD'
SD         -- name of the spherical shell to be plotted
"""

    if not SD.startswith('Sphere'):
        raise NameError('the name of SD={} must start with \'Sphere\''.format(SD))
    a=AdjustGrid[SD] # shortcut

    t=a['Extents'][:,0]
    Nr=a['Extents'][:,1]
    Ntheta=a['Extents'][:,2]

    c1=next(ax._get_lines.prop_cycler)['color']
    c2=next(ax._get_lines.prop_cycler)['color']

    ax.plot(t,(Nr-20.)/10,color=c1, label='(Nr-20)/10')
    ax.plot(t,(Ntheta-20.)/10,color=c2,label='(Ntheta-20)/10')

    tmp=a['Bf0I1']
    ax.plot(tmp[:,0],tmp[:,3],'+-',color=c1,label='TruncErrorExcess - r',linewidth=2)
    tmp=a['Bf1S2']
    ax.plot(tmp[:,0],tmp[:,3],'+-',color=c2,label='TruncErrorExcess - theta',linewidth=2)
    ax.set_xlabel('t/M')
    ax.legend();
    ax.set_title(SD)


def PlotSubdomainConstraints(ax, GhCe, N=5, Ngrey=0):
    """
Make a plot of constraints.
  ax -- axes to plot into
  GhCe -- a dictionary containing constaint info, as obtained with segment_utils.LoadDat_from_segments
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
