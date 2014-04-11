from matplotlib import pyplot as pl
from matplotlib import cm
import numpy as np
from astropy.io import ascii
import sys
import os

thisfile = sys.argv[0]
# when running from ipython, sys.argv[0] points to the python bin dir instead
# of the source directory containing the pipeline code
if 'ipython' in thisfile: thisfile = __file__
thisdir = os.path.dirname(thisfile)

_datfile = 'lensed_galaxies.txt'


def mkTvisHistFigs():

    import plotsetup
    plotsetup.fullpaperfig(1, [8,3])

    datfile1 = 'lensed_galaxies.txt'
    dat1 = ascii.read(datfile1, header_start=-1)
    mu1 = dat1['mu1']
    z1 = dat1['z']
    dt1 = dat1['dt1']

    datfile2 = 'lensed_galaxies2.txt'
    dat2 = ascii.read(datfile2)
    mu2 = dat2['col6']
    z2 = dat2['col5']
    dt2 = dat2['col1']
    dterr2p = dat2['col2']/11.2
    dterr2m = dat2['col3']/11.2
    dterr2 = dat2['col4']/11.2

    ax1 = pl.subplot(131)
    z = np.append( z1, z2 )
    dt = np.append( dt1, dt2 )

    i5 = np.where( (dt>0) & (dt<5) )[0]
    print( len(i5) / float(len(dt) ) )
    print( len(i5) )
    randoff = (np.random.randn( len(z) ) -0.5 ) * 0.1
    z = z+randoff
    pl.hist( z , normed=False, bins=np.arange(0.1,4.51,0.2), color='darkcyan')#, alpha=0.3  )
    # pl.hist( z[i5] , bins=np.arange(0.1,4.51,0.2), color='darkorange')#, alpha=0.3  )
    ax1.set_xlabel('Redshift')
    ax1.set_ylabel(r'Number of Lensed Galaxy Images')

    ax2 = pl.subplot(132)
    pl.hist( dt[i5] , bins=np.arange(0.0,5.01,0.5), color='darkcyan')#, alpha=0.3  )
    ax2.set_xlabel(r'Time Delay $\Delta$t')


    ax3 = pl.subplot(133)

    ipos = np.where( (dt2>0) )
    # err = dterr2[i5] / dt2[i5]
    err = dterr2[ipos] / dt2[ipos]
    pl.hist( err ,normed=False,  bins=np.arange(0.0,3.0,0.05), color='darkcyan' )#, alpha=0.3  )
    ax3.set_xlabel(r'Time Delay \%Error')


    # pl.hist( z[i5] , bins=np.arange(0.1,4.51,0.2), color='darkorange', alpha=0.3  )
    #ax3 = pl.axes( [0.5,0.5,0.4,0.4], transform=ax2.transAxes )
    #ax3.hist( dt, bins=np.arange(0.0,0.1,0.02), color='darkorange', alpha=0.3  )

    fig = pl.gcf()
    fig.subplots_adjust(left=0.08, bottom=0.18,right=0.95, top=0.90, wspace=0.05,)
    
    ax1.set_xlim(0.3,4.6)
    ax2.set_xlim([0.01, 4.99])
    ax3.set_xlim([0.01, 0.99])


    # ax1.set_xticklabels([1,2,3,4])
    ax2.set_yticklabels([])
    ax3.set_yticklabels([])
    # return( i5, dt )

def computeYield( datfile = _datfile ) : 
    """ compute the SN yield per snapshot """
    import cosmo

    tvis1, tvis2 = mkTvisFunctions( etime=20, filter1='f140w', filter2='f110w' )

    # read in the lensed galaxy data
    dat = ascii.read(datfile, header_start=-1)

    cluster = dat['cluster']
    dt = dat['dt1']

    # Limit to galaxies that produce a useful time delay : 
    igood = np.where( (dt>0) & (dt<20) )[0]
    mu = dat['mu1'][igood]
    z = dat['z'][igood]
    dt = dt[igood]
    m850 = dat['m850'][igood]
    m160 = dat['m160'][igood]

    # fix the non-detections at H = 28.5 AB mag
    i99 = np.where( (m160>90) | (m160<19) )
    m160[i99] = 28.5

    i99 = np.where( m850> 90 )
    m850[i99] = 28.5

    NgalTest = len(igood) # total number of galaxies in our 3-cluster test set
    NgalAve = 15 # number of useful galaxy images per snapshot

    tvisIa = np.array( [ tvis1(z1,mu1)[0] for z1,mu1 in zip( z, mu) ] )/365. # years
    tvisCC = np.array( [ tvis2(z2,mu2)[0] for z2,mu2 in zip( z, mu) ] )/365. # years

    Mpc2cm = 3.0857e24  # Mpc to cm conversion

    # convert from observed F160W mag to demagnified B band absolute magnitude
    # (asserting that F160W is approximately rest-frame B band for most of these galaxies)
    muLCDM = cosmo.mu( z ) 
    MB  = m160 + 2.5*np.log10(mu) - muLCDM 
    LBgal = 10**(-0.4*( MB - 5.48 ) - 10 )  # [10^10 Lsun] 
    # and then convert the B band luminosity (LB) to stellar mass M using M/L ratios from Bell+ 2003
    Mgal = 10**(0.5) * LBgal  # [10^10 Msun]

    # Convert the rest-frame UV magnitude into star formation rate
    # First, convert from the measured AB mags to log(flux) and thence to luminosity 
    # (asserting that f850lp is approximately rest-frame UV for most of these galaxies)
    logFuv = (m850 + 2.5*np.log10(mu) + 2.5*np.log10(1+z) + 48.6)/(-2.5) # F in [ergs/s/cm^2/Hz]
    logFB = (m160 + 2.5*np.log10(mu) + 2.5*np.log10(1+z) + 48.6)/(-2.5) # F in [ergs/s/cm^2/Hz]
    Luv =  4 * np.pi * (DL(z)*Mpc2cm)**2 * 10**(logFuv) # [ergs/s/Hz]
    LB =  4 * np.pi * (DL(z)*Mpc2cm)**2 * 10**(logFB) # [ergs/s/Hz]

    # now convert to star formation rate following Kennicutt (1998)
    # and Leitherer (1995)
    # SFR = 1.57 * 1.4e-28 * Luv # [Msun/yr]
    # SFR = 1.57 * 1.4e-28 * LB # [Msun/yr]

    # from dahlen et al 2007, at z=2 using a Salpeter IMF.
    SFR = 9.3e-29 * LB  # [ Msun /yr ] 

    # From Madau+ 1998 and Bouwens+ 2012
    # SFR = 1.25e-28 * LB  # [ Msun /yr ] 

    # now convert from star formation rate to CC SN rate : 
    kcc = 7   # efficiency factor, assuming a Salpeter IMF 
    # consistent with observational constraints, e.g. Dahlen+ 2012
    SNRcc = kcc * 1e-3 * SFR
    
    ## and for the Ia rate use the Scann+Bildsten A+B model: 
    #A = 4.4e-2 # for stellar mass in units of 10^10 Msun
    #B = 0.26  # for SFR in units of Msun/yr
    #SNRIa = (A * Mgal + B * SFR ) / 100. # SNIa / yr

    # use the measured volumetric SNIa rate at z~2 from Rodney+ 2014
    # and convert to SNuB using the measured B band luminosity density
    # at z~2 from Reddy&Steidel 2009
    SNRIa_vol = 0.5e-4 # SNIa / yr / Mpc^3  (R14)
    rho_LB = 1.4e27  # ergs/s/Hz/Mpc^3  in UV (R&S 2009)
    Lsun  = 1.9e33   # ergs/s
    nuB = 6.81e14 # Hz
    # SNRIa = (SNRIa_vol / rho_SFR) * Luv  # SNIa / yr

    rho_Lsun = rho_LB*nuB / (1e10*Lsun)  # 10^10 Lsun/Mpc^3
    SNRIa = (SNRIa_vol / (rho_Lsun))* LBgal  # SNIa / yr

    NsnIa = ( SNRIa * tvisIa / (1+z) ).sum() / NgalTest * NgalAve
    NsnCC = ( SNRcc * tvisCC / (1+z) ).sum() / NgalTest * NgalAve

    print( "Expected SN yield per cluster snapshot")
    print( "  computed using observed volumetric high-z SN rates for Ia and SFR scaling for CC" )
    print( "NsnIa  = %.3e"%( NsnIa ) )
    print( "NsnCC  = %.3e"%( NsnCC ) )

    print( "==> Expected total yield per 100 snapshots: ")
    print( "Nsn  = %.3f"%( (NsnIa+NsnCC)*100. ) )

    # alternative approach, using rest-frame B-band luminosity
    # and the SNuB prescriptions from Li+ 2011 for an Scd galaxy
    SNRuBIa =  (0.2/100) * (LBgal/2)**(-0.25)  # SN / yr / 10^10 Lsun
    SNRuBcc =  (0.7/100) * (LBgal/2)**(-0.25)  # SN / yr / 10^10 Lsun

    # adopting the SNR per unit mass from Mannucci+ 2005
    # for Irregular galaxies
    # SNRuBIa = 0.8 / 100 # SN /  yr / 10^10 Msun  
    # SNRuBcc = 2.25 / 100 # SN /  yr / 10^10 Msun 

    # SNRuBIa = 0.2 / 100 # SN /  yr / 10^10 Msun  
    # SNRuBcc = 0.8 / 100 # SN /  yr / 10^10 Msun 

    SNRIa = SNRuBIa * LBgal 
    SNRcc = SNRuBcc * LBgal 

    NsnIa = ( SNRIa * tvisIa / (1+z) ).sum() / NgalTest * NgalAve
    NsnCC = ( SNRcc * tvisCC / (1+z) ).sum() / NgalTest * NgalAve
    
    print( "computed using Li+ 2011 SNuB and galaxy rest-frame B band luminosity" )
    print( "NsnIa  = %.3e"%( NsnIa ) )
    print( "NsnCC  = %.3e"%( NsnCC ) )

    print( "==> Expected total yield per 100 snapshots: ")
    print( "Nsn  = %.3f"%( (NsnIa+NsnCC)*100. ) )
    
    

    

def mkTvisFunctions( etime=20, filter1='f140w', filter2='f110w' ):
    """Read in simulations, build 2-D interpolators and return them as
    two functions that give the expected mean visibility time for a Ia
    and a CC at a given  (z,mu) position.
    """
    import os
    from scipy import interpolate as scint
    from astropy.io import ascii

    tvisdatfile1 = 'snIa_%s_%02imin_tvis.dat'%(filter1,etime)
    tvisdatfile2 = 'snII_%s_%02imin_tvis.dat'%(filter2,etime)
    tvisdat2 = ascii.read( tvisdatfile2, header_start=-1 )

    interpolators = []
    for tvisdatfile in [tvisdatfile1,tvisdatfile2] : 
        tvisdat = ascii.read( tvisdatfile, header_start=-1 )
        zgrid = tvisdat['z']
        mugrid =  np.array([ int(col[-2:]) for col in tvisdat.colnames
                             if col.startswith('tvis') ] )
        tvisgrid = np.array( [ tvisdat[col] for col in tvisdat.colnames
                               if col.startswith('tvis') ] )
        tvisinterp = scint.interp2d( zgrid, mugrid, tvisgrid, bounds_error=False,
                                     fill_value=None ) # extrapolate!
        interpolators.append( tvisinterp )
    return( interpolators )


def DL( z ) :
    """ approximate luminosity distance in Mpc for a flat "737" LCDM universe 
    Using the approximation of Wickramasinghe and Ukwatta (2010, arXiv:1003.0483)
    """
    H0=70
    Om=0.3
    Ol=0.7

    c = 299792.458
    alpha = 1 + 2*Ol/(1-Ol)/(1+z)**3
    x = np.log( alpha + np.sqrt( alpha*alpha - 1 ) )

    x0 = 2.42 
    psi0 = -2.210
    psi = lambda y : 3 * y**(1/3.) * 2**(2/3.) * ( 1 - (y**2/252.) + (y**4/21060) ) + psi0
    
    dL = c / (3*H0) * (1+z) / (Ol**(1/6.)*(1-Ol)**(1/3.)) * ( psi(x0) - psi(x) )
    return( dL )


def computeLsun( datfile = _datfile ):
    import cosmo


    dat = ascii.read(datfile, header_start=-1)

    cluster = dat['cluster']
    dt = dat['dt1']
    mu = dat['mu1']
    z = dat['z']
    m850 = dat['m850']
    m160 = dat['m160']

    # fix the non-detections at H = 28.5 AB mag
    i99 = np.where( m160< 99 )
    m160[i99] = 28.5

    # Limit to galaxies that produce a useful time delay : 
    igood = np.where( (dt>0) & (dt<20) )[0]

    dl = cosmo.DLFw( z ) # luminosity distance, assuming a flat universe
    dm_mu = 2.5*np.log10( mu ) # lensing magnification in magnitudes

    MB = m160 - dl + dm_mu # de-magnified F160W mag in AB mags /arcsec2 
    #(not actually B band, but since most are at z~2, this is not too far off, maybe?)

    # MBsun = 4.74 + 0.65 # abs mag of sun 
    # MB = 27.05 - 2.5 log10( I ) 
    # Lsun = 3.839e33 # erg/s

    LBsun10 = (-0.4*(MB-27.05) ) # gal. luminosity in units of 10^10 Lsun/arcsec2

    # at z~2, 1" ~ 5 kpc, so lets say our galaxies are ~10 kpc on a side so ~2x2"
    # sanity check: a typical unlensed galaxy at z~2 is like the primo or wilson host,
    #  so it subtends about 10x10 pixels at 0.09"/pix ~ 1 arcsec^2
    # So it seems reasonable to make ~no correction for angular size here. 
    pl.hist( LBsun10, color='b', alpha=0.3, bins=np.arange(-1,7,0.5), label='All galaxy images' )
    pl.hist( LBsun10[igood], color='r', alpha=0.3, bins=np.arange(-1,7,0.5), label='with time delay to next image < 20 yrs ')

    pl.xlabel('Lensed galaxy luminosity [$10^{10} L_{\odot}$]')
    pl.ylabel('Number of lensed galaxy images')
    pl.legend(loc='upper right')
    return( LBsun10[igood] )



def maghist(datfile=_datfile, expand=False, debug=False):
    if debug:
        import pdb
        pdb.set_trace()

    # Read in the A1689 multiply imaged galaxy catalog
    dat = ascii.read(datfile)
    igal = dat['gal']
    iim = dat['im']
    z = np.abs(dat['z'])
    izphot = np.where(dat['z']<0)[0]
    izspec = np.where(dat['z']>0)[0]
    mu = dat['mu']
    itoohigh = np.where(mu>25)[0]
    mu[itoohigh] = np.random.random(len(itoohigh))*10+20
    idx = dat['idx'] - 1

    if expand :
        # Expand the A1689 catalog by a factor of 4, applying random offsets to
        # get a mock catalog representing the full sample
        zoffsets = np.random.normal(loc=0.0, scale=0.3, size=len(z)*4)
        z4 = np.array( 4 * z.tolist() ) + zoffsets
        muoffsets = np.random.normal(loc=0.0, scale=1.0, size=len(mu)*4)
        mu4 = np.abs( np.array(4 * mu.tolist() ) + muoffsets )
        idx4 = np.array( 4 * idx.tolist() )
        igal4 = np.array( 4 * igal.tolist() )
        z = z4
        mu = mu4
        idx = idx4
        igal = igal4

    imaxmags = []
    iminmags = []
    for ig in np.unique(igal):
        idx_thisgal = idx[np.where(igal == ig)]
        mu_thisgal = mu[idx_thisgal]
        imaxmags.append(idx_thisgal[np.argmax(mu_thisgal)])
        iminmags.append(idx_thisgal[np.argmin(mu_thisgal)])

    imaxmag = np.array( [ np.argmax(idx[np.where(igal == igal)])
                          for igal in np.unique(igal) ] )
    iminmag = np.array( [ np.argmin(mu[np.where(igal == igal)])
                        for igal in np.unique(igal)])
    # mu12diff = mu[imaxmags] - mu[iminmags]

    pl.clf()

    pl.subplot(1,3,1)
    pl.plot(z[iminmags], mu[iminmags], color='b', marker='o', alpha=0.3,
            ms=10, ls=' ')
    pl.plot(z[imaxmags], mu[imaxmags], color='r', marker='s', alpha=0.3,
            ms=10, ls=' ')

    pl.subplot( 1, 3, 2 )
    maxmaghist2d, xedges, yedges = np.histogram2d(z[imaxmags], mu[imaxmags],
                                                  bins=15)
    minmaghist2d, xedges, yedges= np.histogram2d(z[iminmags], mu[iminmags],
                                                 bins=15)
    # pl.contour(z[imaxmags], mu[imaxmags], maxmaghist2d )
    pl.contour( xedges[:-1], yedges[:-1],
                maxmaghist2d, colors='r',
    )
    # 4,2,0] )
    # pl.contour( minmaghist2d, colors='b', levels=[4,2,0] )

    # pl.hist2d(z[iminmags], mu[iminmags], bins=15)


    pl.subplot(1,3,3)
    pl.hist2d(z[imaxmags], mu[imaxmags], bins=15)

def plot_tvis_lines( snIadatfile='snIa_tvis.dat', snIIdatfile='snII_tvis.dat'  ):
    """ Plot SN visibility time vs redshift for a range of mu values using
    .dat tables generated by snapsim.py
    """
    dat1 = ascii.read( snIadatfile, header_start=-1 )
    dat2 = ascii.read( snIIdatfile, header_start=-1 )

    mulist = [2,4,6,10,15,20]
    mucolorlist = ['m','b','c','g','r','k']

    pl.clf()
    ax1 = pl.subplot( 1,2,1 )
    ax2 = pl.subplot( 1,2,2, sharex=ax1, sharey=ax1 )

    for dat,ax,sntype in zip( [dat1,dat2], [ax1,ax2], ['Ia','II']) :

        z = np.array( dat['z'] )
        for mu, mucolor in zip(mulist,mucolorlist) :
            tvis = np.array( dat['tvis%02i'%mu] )
            err = np.array( dat['err%02i'%mu] )
            tvismax = tvis + err
            tvismin = np.max( [np.zeros(len(tvis)), tvis-err], axis=0 )
            # ax.fill_between( z, tvismin, tvismax, color=mucolor, alpha=0.3 )
            ax.plot( z, tvis, marker=' ', ls='-', color=mucolor, label='%i'%mu )
            z10 = z[ np.where(tvis<12)[0][0] ]
            ax.text( z10,10, '%i'%mu, color=mucolor, ha='center', va='center',
                     backgroundcolor='w' )

        # ax.legend(loc='upper right')
        ax.set_xlabel('Redshift')
        ax.set_ylabel('Visibility Time [days]')
        ax.text(0.95,0.95,'Type %s SN'%sntype, ha='right',va='top',
                transform=ax.transAxes, fontsize='large' )
        ax.set_ylim( 0, 140 )
        ax.set_xlim( 0.8, 3.2 )
        ax.text(1.0,10,'$\mu$=',ha='right',va='center',
                backgroundcolor='w' )

    fig = pl.gcf()
    fig.subplots_adjust( left=0.12, right=0.88, bottom=0.12, top=0.95, wspace=0 )

    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.set_ticks_position('both')
    ax2.set_ylabel('Visibility Time [years]', rotation=-90 )

    # ax1.set_xlim(0.9,3.2)

    ax2.set_yticks( np.array([0.1,0.2,0.3])*365 )

    ax2.set_yticklabels( [0.1,0.2,0.3] )

    ax1.set_ylim(0,120)

    return( dat1, dat2 )


def mk_tvis_fig( obsdat=_datfile ):
    import plotsetup


    plotsetup.halfpaperfig( figsize=[5,5] )
    ax1 = pl.subplot( 111 )

    # dat1 = 'snIa_F140W_12min_tvis.dat'
    # plot_tvis_contours( dat=dat1, obsdat=obsdat, points=points, ls='dashdot' )

    dat1 = 'snIa_F140W_20min_tvis.dat'
    plot_tvis_contours( dat=dat1, obsdat=obsdat, points=True, ls='solid' )

    # dat1 = 'snII_F110W_20min_tvis.dat'
    # plot_tvis_contours( dat=dat1, obsdat=obsdat, points=True, ls='dashed' )

    # dat1 = 'snIa_F140W_30min_tvis.dat'
    # plot_tvis_contours( dat=dat1, obsdat=obsdat, points=False, ls='solid' )
    # plot_tvis_contours( dat=dat3, obsdat=obsdat, expand=expand, points=points, ls='dashed' )
    ax1.text( 0.05, 0.95, 'Type Ia\n Visibility\n Windows\n for 20 min \n F140W snaps', transform=ax1.transAxes, ha='left', va='top', fontsize='large' )
    #ax1.text( 1.85,22.4,r'130 days', rotation=75, color='darkred', backgroundcolor='w', ha='center', fontsize=12 )
    ax1.text( 2.9,25.,'t$_{vis}$=100 days', va='bottom', ha='right', rotation=0, color='darkred', fontsize=12 )
    ax1.text( 3.1,25.,'50', va='bottom', ha='center', rotation=0, color='darkcyan', fontsize=12 )
    ax1.text( 3.32,25.,'10', va='bottom', ha='center', rotation=0, color='darkorchid', fontsize=12 )
    # ax1.text( 2.97,22.4,'40', ha='center', rotation=85, color='darkblue', backgroundcolor='w',fontsize=12 )


    # # ax2 = pl.subplot( 122 )
    # # plot_tvis_contours( dat=dat2, obsdat=obsdat, expand=expand, points=points, ls='solid' )
    # plot_tvis_contours( dat=dat2, obsdat=obsdat, points=points, ls='solid' )
    # # ax2.text( 0.05, 0.95, 'Type II', transform=ax2.transAxes, ha='left', va='top', fontsize='large' )
    # ax2.text( 0.05, 0.95, 'Type II\n Visibility\n Windows\n for %02i min \n F110W snaps'%etime, transform=ax2.transAxes, ha='left', va='top', fontsize='large' )
    # ax2.text( 1.72,22.4,'130 days', rotation=74, color='darkred', backgroundcolor='w', ha='center', fontsize=12 )
    # ax2.text( 1.87,22.4,'100', ha='center', rotation=72, color='darkorange', backgroundcolor='w',fontsize=12 )
    # ax2.text( 2.25,22.4,'70', ha='center', rotation=74, color='darkcyan', backgroundcolor='w',fontsize=12 )
    # ax2.text( 2.75,22.4,'40', ha='center', rotation=76, color='darkblue', backgroundcolor='w',fontsize=12 )
    # ax2.text( 3.1, 22.4,'10', ha='center', rotation=78, color='darkorchid', backgroundcolor='w',fontsize=12 )
    #
    # ax2.yaxis.set_ticks_position('right')
    # ax2.yaxis.set_ticks_position('both')
    # ax2.yaxis.set_label_position('right')
    # ax2.set_ylabel( 'Magnification, $\mu$', labelpad=20, rotation=-90 )
    fig = pl.gcf()
    fig.subplots_adjust( left=0.12,bottom=0.13, right=0.92, top=0.9, wspace=0.05 )
    pl.draw()

def plot_tvis_contours( dat='snIa_tvis.dat', obsdat=_datfile, 
                        points=True, ls='solid' ):
    """ Plot SN visibility time contours in the mu vs redshift plane
    """

    if isinstance( dat, str ) :
        dat = ascii.read( dat, header_start=-1 )

    # Read in the grid of simulated visibility times
    zgrid = dat['z']
    mugrid =  np.array([ int(col[-2:]) for col in dat.colnames
                     if col.startswith('tvis') ] )
    tvismatrix = np.array( [dat['tvis%02i'%mu] for mu in mugrid ] )


    if points :
        # Read in the observed multiply imaged galaxy catalog
        if isinstance( obsdat, str ) :
            obsdat = ascii.read(obsdat, header_start=-1)
        zobs = np.abs(obsdat['z'])  # negative z's are photoz
        muobs = obsdat['mu1']
        dt = obsdat['dt1']
        tvisIa = obsdat['tvisIa']
        tvisII = obsdat['tvisII']

        # isolate the images that can provide a useful time delay
        #   - time delay to next image < 20 years
        #   - visibility time > 0
        # igood = np.where( (dt>0) & (dt<20) & ((tvisIa>0) | (tvisII>0)) )[0]
        igood = np.where( (dt>0) & (dt<20) )[0]
        muTD = np.min( [muobs[ igood ], 20 + np.random.random_sample(len(igood))*5], axis=0 )
        zTD = zobs[ igood ]

        pl.plot(zTD, muTD, color='k', marker='o', alpha=0.3,
                ms=10, ls=' ')

    ax = pl.gca()
    ax.set_xlim( 0.95, 3.6 )
    ax.set_ylim( 1.01, 25 )

    colorlist = ['darkorchid','darkcyan','darkred']
    levels = [10,50,100]
    CS = pl.contour( zgrid, mugrid, tvismatrix, levels=levels, colors=colorlist, linestyles=ls )

    ax.set_xlabel('Source Redshift')
    ax.set_ylabel('Magnification, $\mu$')

    # pl.clabel(CS, levels, inline=True, fontsize=14, fmt='%i' )
    # return( zmuhistMI, xedges, yedges )




def print_tvis_means( tvisfile='lensed_galaxies_tvis.dat' ):
    """Read in the tvis data for all the lensed galaxies
    in our observed data set. Compute and print the mean
    visibility time for Ia and II for 12,20,30-min snaps.
    """
    from astropy.io import ascii
    tvisdat = ascii.read( tvisfile )
    pl.clf()

    # TODO : isolate the time-delay capable images

    igal = tvisdat['gal']
    idx = tvisdat['idx'] - 1
    muobs = tvisdat['mu']

    iTD = []
    iLast = []
    for ig in np.unique(igal):
        idx_thisgal = idx[np.where(igal == ig)]
        ithisgalLast = np.argsort( muobs[idx_thisgal] )
        iLast.append( idx_thisgal[ithisgalLast[-1]])
        iTD += idx_thisgal[ithisgalLast[:-1]].tolist()

    # hist12, binedges = np.histogram( tvisIa_12, bins=np.arange(1,150,10) )
    # hist20, binedges = np.histogram( tvisIa_20, bins=np.arange(1,150,10) )
    # hist30, binedges = np.histogram( tvisIa_30, bins=np.arange(1,150,10) )
    # pl.bar( binedges[:-1], hist12, width=10, color='b', alpha=0.3, linewidth=1 )
    # pl.bar( binedges[:-1], hist20, width=10, color='g', alpha=0.3, linewidth=1 )
    # pl.bar( binedges[:-1], hist30, width=10, color='r', alpha=0.3, linewidth=1 )

    # pl.plot( binedges[:-1], hist12, drawstyle='steps-mid', color='b', linewidth=1 )
    # pl.plot( binedges[:-1], hist20, drawstyle='steps-mid', color='g', linewidth=1 )
    # pl.plot( binedges[:-1], hist30, drawstyle='steps-mid', color='r', linewidth=1 )

    print( 'ACS/50: %02i days  %02i days'%( np.mean( tvisdat['tvisIa_50'][iTD] ), np.mean( tvisdat['tvisII_50'][iTD] ) ) )
    print( '45min : %02i days  %02i days'%( np.mean( tvisdat['tvisIa_45'][iTD] ), np.mean( tvisdat['tvisII_45'][iTD] ) ) )
    print( '30min : %02i days  %02i days'%( np.mean( tvisdat['tvisIa_30'][iTD] ), np.mean( tvisdat['tvisII_30'][iTD] ) ) )
    print( '20min : %02i days  %02i days'%( np.mean( tvisdat['tvisIa_20'][iTD] ), np.mean( tvisdat['tvisII_20'][iTD] ) ) )
    print( '12min : %02i days  %02i days'%( np.mean( tvisdat['tvisIa_12'][iTD] ), np.mean( tvisdat['tvisII_12'][iTD] ) ) )


    # return( tvisIa_12, tvisIa_20, tvisIa_30 )





def mk_lightcurve_fig(  ):
    from hstsnpipe.tools import snana
    from hstsnpipe.tools.figs import plotsetup
    plotsetup.fullpaperfig( figsize=[8,4] )
    pl.clf()

    simIa = snana.SimTable( 'snIa_zgrid' )
    FLTMATRIX = simIa.FLT.reshape( simIa.LCMATRIX.shape )
    bandlist = FLTMATRIX[0,0,0,0,:,0]
    def iband(band) : 
        return( np.where( bandlist==band )[0][0] )
    iz18 = np.argmin( np.abs( simIa.z-1.8 ) )
    ilp0 = np.argmin( np.abs( simIa.x1 ) )
    iF160W = iband('H')
    iF140W = iband('N')
    iF110W = iband('M')
    muIa = 5.

    ax1 = pl.subplot(121)
    ax2 = pl.subplot(122, sharex=ax1 )
    # snIaMAGH = simIa.LCMATRIX[ ilp0, 0, 0, iz18, iF160W, : ]
    snIaMAGN = vega2ab( simIa.LCMATRIX[ ilp0, 0, 0, iz18, iF140W, : ], 'F140W' )
    snIaMAGM = vega2ab( simIa.LCMATRIX[ ilp0, 0, 0, iz18, iF110W, : ], 'F140W' )
    snIaMJD = simIa.TOBS[ iz18, : ]

    # ax1.plot( snIaMJD, snIaMAGH - 2.5*np.log10( muIa ), 'k-', lw=2, marker=' ', )
    ax1.plot( snIaMJD, snIaMAGN - 2.5*np.log10( muIa ), 'r-', lw=2, marker=' ', )
    ax2.plot( snIaMJD, snIaMAGM - 2.5*np.log10( muIa ), 'r-', lw=2, marker=' ', )
    ax1.invert_yaxis()

    simII = snana.SimTable( 'snII_zgrid' )
    FLTMATRIX = simII.FLT.reshape( simII.LCMATRIX.shape )
    bandlist = FLTMATRIX[0,0,0,0,:,0]
    def iband(band) : 
        return( np.where( bandlist==band )[0][0] )
    iz20 = np.argmin( np.abs( simII.z-2.0 ) )
    iF140W = iband('N')
    iF110W = iband('M')
    muII = 12.
    ilpII=13  # select one representative II-P model

    snIIMAGN = vega2ab( simII.LCMATRIX[ ilpII, 0, 0, iz20, iF140W, : ], 'F140W' )
    snIIMAGM = vega2ab( simII.LCMATRIX[ ilpII, 0, 0, iz20, iF110W, : ], 'F110W' )
    snIIMJD = simII.TOBS[ iz20, : ]
    ax1.plot( snIIMJD, snIIMAGN - 2.5*np.log10( muII ), 'b--', lw=2, marker=' ', )
    ax2.plot( snIIMJD, snIIMAGM - 2.5*np.log10( muII ), 'b--', lw=2, marker=' ', )

    ax1.axhline( vega2ab( 24.6, 'f140w'), color='k',ls=':', lw=2, )
    ax2.axhline( vega2ab( 25.3, 'f110w'), color='k',ls=':', lw=2, )

    # ax2.text(0.95,0.95,'F160W', fontsize='large', ha='right',va='top', transform=ax2.transAxes )
    ax1.text(0.95,0.95,'F140W', fontsize='large', ha='right',va='top', transform=ax1.transAxes )
    ax2.text(0.95,0.95,'F110W', fontsize='large', ha='right',va='top', transform=ax2.transAxes )

    ax1.text( 115, 25.65, '30 min SNAP\n detection limit', ha='right',va='bottom')
    # ax2.text( 115, 26.0, '30 min SNAP\n detection limit', ha='right',va='bottom')

    ax1.text( 38, 25.3, 'Type Ia SN\n z=1.8, $\mu$=5\n t$_{vis}$=75 days', ha='left',va='bottom', color='r', fontsize='large')
    # ax2.text( 25, 25.5, 'Type II-P SN\n z=2.0, $\mu$=12\n t$_{vis}$=50 days', ha='left',va='center', color='b', fontsize='large')
    # ax1.text( 67, 25.6, 'Type II-P SN\n z=2.0, $\mu$=12\n t$_{vis}$=50 days', ha='left',va='top', color='b', fontsize='large', backgroundcolor='w' )
    ax2.text( 67.5, 26.15, 'Type II-P SN\n z=2.0, $\mu$=12\n t$_{vis}$=50 days', ha='left',va='top', color='b', fontsize='large', backgroundcolor='w' )

    ax1.set_xlim(-35,119)
    ax1.set_ylim(27, 27-2.51)
    ax2.set_ylim(27.4, 27.4-2.51)
    ax2.yaxis.set_ticks_position('right')
    ax2.yaxis.set_ticks_position('both')
    ax2.yaxis.set_label_position('right')
    fig = pl.gcf()
    fig.subplots_adjust( left=0.13, right=0.87, bottom=0.12, top=0.95, wspace=0.05 )
    ax2.set_xlabel('Observer-frame Time [days]')
    ax1.set_xlabel('Observer-frame Time [days]')
    ax1.set_ylabel('Observed F140W magnitue [AB]')
    ax2.set_ylabel('Observed F110W magnitue [AB]', labelpad=20, rotation=-90)


ZPT_WFC3_IR_AB = {'F105W':26.0974,
                  'F110W':26.6424,
                  'F125W':26.0449,
                  'F140W':26.2608,
                  'F160W':25.7551,
                  'F098M':25.5041,
                  'F127M':24.4545,
                  'F139M':24.2880,
                  'F153M':24.2725 }

ZPT_WFC3_IR_VEGA = {'F105W':25.4523,
                    'F110W':25.8829,
                    'F125W':25.1439,
                    'F140W':25.1845,
                    'F160W':24.5037,
                    'F098M':24.9424,
                    'F127M':23.4932,
                    'F139M':23.2093,
                    'F153M':23.0188 }

def vega2ab( magVega, filter ):
    zptAB = ZPT_WFC3_IR_AB[filter.upper()]
    zptVega = ZPT_WFC3_IR_VEGA[filter.upper()]
    return( magVega  - zptVega + zptAB )


def plot_pkmags_AB( simIa ):

    Hpk = vega2ab( simIa.LCMATRIX[7,0,0,:,4,19], 'F160W' )
    Npk = vega2ab( simIa.LCMATRIX[7,0,0,:,3,19], 'F140W' )
    Jpk = vega2ab( simIa.LCMATRIX[7,0,0,:,2,19], 'F125W' )
    Mpk = vega2ab( simIa.LCMATRIX[7,0,0,:,1,19], 'F110W' )
    Ypk = vega2ab( simIa.LCMATRIX[7,0,0,:,0,19], 'F105W' )

    pl.clf()
    ax = pl.gca()

    ax.plot( simIa.z, Hpk, 'r-' )
    ax.plot( simIa.z, Npk, color='darkorange', ls='-', marker=' ')
    ax.plot( simIa.z, Jpk, 'g-' )
    ax.plot( simIa.z, Mpk, 'b-' )
    ax.plot( simIa.z, Jpk, 'g-' )
    ax.plot( simIa.z, Mpk, 'b-' )
    ax.plot( simIa.z, Ypk, color='darkorchid', ls='-', marker=' ' )

    ax.set_xlabel('Redshift')
    ax.set_ylabel('SN Ia peak magnitude (AB)')
    ax.invert_yaxis()
    ax.set_xlim(0.9,3.49)
    ax.set_ylim(32.5,24.5)

    Htxt = ax.text( 3.3, 28.9, 'F160W', backgroundcolor='w', ha='center', va='center', color='r' )
    Ntxt = ax.text( 3.3, 29.6, 'F140W', backgroundcolor='w', ha='center', va='center', color='darkorange' )
    Jtxt = ax.text( 3.3, 31.4, 'F125W', backgroundcolor='w', ha='center', va='center', color='g' )
    Mtxt = ax.text( 3.3, 32.0, 'F110W', backgroundcolor='w', ha='center', va='center', color='b' )
    Mtxt = ax.text( 2.9, 32.2, 'F105W', backgroundcolor='w', ha='center', va='center', color='darkorchid' )


def print_Nleading( obsdat=_datfile):
    """Print the number of leading images from our catalog of
    multiply-imaged galaxies. i.e. count up all the instances of
    multiple images in each system, except the last one, from which we
    could not measure a time delay.
    """
    # Read in the A1689 multiply imaged galaxy catalog
    if isinstance( obsdat, str ) :
        obsdat = ascii.read(obsdat)
    igal = obsdat['gal']
    zobs = np.abs(obsdat['z'])  # negative z's are photoz
    muobs = obsdat['mu']
    idx = obsdat['idx'] - 1

    # isolate the leading images (any images that
    # have later images from which a time delay
    # could be measured.
    muTD, zTD = [], []
    muLast, zLast = [], []
    for ig in np.unique(igal):
        idx_thisgal = idx[np.where(igal == ig)]
        mu_thisgal = sorted( muobs[idx_thisgal] )
        muTD += mu_thisgal[:len(mu_thisgal)-1]
        zTD += zobs[idx_thisgal][:len(mu_thisgal)-1].tolist()
        muLast.append( mu_thisgal[-1] )
        zLast.append( zobs[idx_thisgal][-1] )
    Ntot = len(muobs)
    Nlead = len(muTD)
    print( '%i leading images out of %i total images (%.1f pct)'%(
        Nlead, Ntot, 100.*Nlead/Ntot ) )


def mk_dt_hist_fig( datfile=_datfile ):
    """ Make a figure showing the histogram of time delays
    for our example set, limiting to just those lensed galaxy
    images that are suitable for SN time delay measurements
    i.e.  mu and z predict tvis > 0 for 30min snap
          not the last image
          time delay to next image < 20 years
    """
    from astropy.io import ascii
    from matplotlib import pyplot as pl

    dat = ascii.read( datfile, header_start=-1 )
    cluster = dat['cluster']
    dt1 = dat['dt1']
    dt2 = dat['dt2']

    pl.clf()
    pl.figure(1)
    pl.hist( dt1, bins=np.arange(0,21,1), color='blue', alpha=0.3, label='LTM' )
    pl.hist( dt2, bins=np.arange(0,21,1), color='red', alpha=0.3, label='NFW')
    pl.ylabel('number of lensed galaxy images')
    pl.xlabel('time delay to next image [years]')
    pl.legend(loc='upper right')

    return( cluster, dt1, dt2 )





def main(argv):
    if argv is None :
        argv == sys.argv

    maghist()

if __name__ == "__main__" :
    main(None)
