from hstsnpipe.tools import snana
import numpy as np
from matplotlib import pyplot as pl

def do_grid_sims_Ia( clobber=False ) :
    """ Execute the SNANA SN Ia simulation : 
    26 steps in x1 (to match the Type II template set)
    50 steps in z 
    50 steps in trest
    """
    snana.simulate.mkinputGrid('snIa_zgrid', inputfile='snIa_zgrid.input',
                               simlibfile='SNAPzgrid.simlib',
                               survey='HST', simtype='Ia',
                               GENFILTERS='YMJNH',
                               ngrid_trest=60, genrange_trest=[-19, 40],
                               ngrid_logz=40, genrange_redshift=[0.9, 3.5],
                               ngrid_lumipar=15, genrange_salt2x1=[-2.5, 2.5],
                               ngrid_colorpar=1, genrange_salt2c=[0, 0],
                               ngrid_colorlaw=1, genrange_rv=[3.1, 3.1],
                               clobber=clobber )
    snana.simulate.mkgridsimlib('SNAPzgrid.simlib', survey='HST',
                                field='default', bands='YMJNH',
                                clobber=clobber)
    snana.simulate.dosim('snIa_zgrid.input', perfect=True )
    simtable = snana.SimTable( 'snIa_zgrid' )
    return( simtable )

def do_grid_sims_II( clobber=False ) :
    """ Execute the SNANA SN II simulation (II-P and IIn) :
    26 steps in lumipar (defined by the Type II template set)
    50 steps in z 
    50 steps in trest
    """
    snana.simulate.mkinputGrid('snII_zgrid', inputfile='snII_zgrid.input',
                               simlibfile='SNAPzgrid.simlib',
                               survey='HST', simtype='II',
                               GENFILTERS='YMJNH',
                               ngrid_trest=60, genrange_trest=[-19, 40],
                               ngrid_logz=40, genrange_redshift=[0.9, 3.5],
                               ngrid_lumipar=1, genrange_salt2x1=[0, 0],
                               ngrid_colorpar=1, GENRANGE_AV=[0,0],
                               ngrid_colorlaw=1, genrange_rv=[3.1, 3.1],
                               clobber=clobber )
    snana.simulate.mkgridsimlib('SNAPzgrid.simlib', survey='HST',
                                field='default', bands='YMJNH',
                                clobber=clobber)
    snana.simulate.dosim('snII_zgrid.input', perfect=True )
    simtable = snana.SimTable( 'snII_zgrid' )
    return( simtable )

def do_single_sims():
    snana.simulate.doSingleSim(simname='snIaTest', z=1, pkmjd=55500, model='Ia', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='JHW', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIa10 =  snana.simulate.doSingleSim(simname='snIa10', z=1.0, pkmjd=55500, model='Ia', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIa15 =  snana.simulate.doSingleSim(simname='snIa15', z=1.5, pkmjd=55500, model='Ia', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIa20 =  snana.simulate.doSingleSim(simname='snIa20', z=2.0, pkmjd=55500, model='Ia', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIa30 =  snana.simulate.doSingleSim(simname='snIa30', z=3.0, pkmjd=55500, model='Ia', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIa50 =  snana.simulate.doSingleSim(simname='snIa50', z=5.0, pkmjd=55500, model='Ia', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)

    snIIP10 =  snana.simulate.doSingleSim(simname='snIIP10', z=1.0, pkmjd=55500, model='201', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIIP15 =  snana.simulate.doSingleSim(simname='snIIP15', z=1.5, pkmjd=55500, model='201', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIIP20 =  snana.simulate.doSingleSim(simname='snIIP20', z=2.0, pkmjd=55500, model='201', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIIP30 =  snana.simulate.doSingleSim(simname='snIIP30', z=3.0, pkmjd=55500, model='201', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)
    snIIP50 =  snana.simulate.doSingleSim(simname='snIIP50', z=5.0, pkmjd=55500, model='201', trestrange=[-15, 30], Av=0, mB=0, x1=0, c=0, dustmodel='mid', survey='HST', field='default', bands='MN', Nobs='max', cadence='min', mjdlist=[], bandlist=[], perfect=True, verbose=False, clobber=False, debug=False)

def mk_tvis_tables_IaII( clobber=False ) :
    """ Run the simulations, read in the simulations, and make the tables of
    visibility time vs redshift vs mu.
    """
    simIa =  do_grid_sims_Ia( clobber=clobber )
    simII = do_grid_sims_II( clobber=clobber )

    for etime in [12,20,30] :
        mk_tvis_table( simIa, filter='F140W', etime=etime, showplots=False )
        mk_tvis_table( simII, filter='F110W', etime=etime, showplots=False )


def update_obstable( infile='lensed_galaxies.txt', outfile='lensed_galaxies_tvis.dat' ):
    """ Read in the mu and z for a set of observed multiply-imaged galaxies from
    the infile.  For each galaxy image, compute the visibility time of a SN
    in that galaxy for 12, 20 and 30-minute snapshots.  Record the results
    to outfile.
    """
    from astropy.io import ascii
    from astropy.table import Column

    from scipy import interpolate as scint

    indat = ascii.read( infile )
    zobs = np.abs( indat['z'] )
    muobs = indat['mu']

    # for each of the 12, 20, and 30-minute snapshots,
    # construct a 2d interpolator for estimating tvis in the mu-z plane
    for sntype, filter in ['Ia','F140W'],['II','F110W'] :
        for etime in [12,20,30] :
            tvisdatfile = 'sn%s_%s_%02imin_tvis.dat'%(sntype,filter,etime)
            tvisdat = ascii.read( tvisdatfile, header_start=2 )
            zgrid = tvisdat['z']
            mugrid =  np.array([ int(col[-2:]) for col in tvisdat.colnames
                                 if col.startswith('tvis') ] )
            tvisgrid = np.array( [ tvisdat[col] for col in tvisdat.colnames
                                   if col.startswith('tvis') ] )
            tvisinterp = scint.interp2d( zgrid, mugrid, tvisgrid, bounds_error=False,
                                         fill_value=None ) # extrapolate!


            tvisobs = np.array( [ tvisinterp(z,mu)[0] for z,mu in zip( zobs, muobs) ] )
            tviscol = Column( tvisobs, name='tvis%s_%02i'%(sntype,etime) )
            indat.add_column( tviscol )

    indat.write(outfile,format='ascii.commented_header')
    print( "Updated table of lensed galaxies with t_vis printed to %s"%outfile)


def mk_tvis_table( sim, filter='F140W', etime=30, showplots=True, verbose=False ) :
    """Create a .dat file giving the visibility time as a function of
    redshift and magnification, for the given filter and exposure time
    (in minutes), using the simulated SN light curves in the given
    simulation table 'sim'.
    """
    import pdb
    band = filter2band( filter )
    maglimit = snapshot_threshold( filter, etime )

    outfile = sim.simname.split('_')[0] + '_%s_%02imin_tvis.dat'%(filter,etime)
    z = sim.z
    # iband = sim.BANDS.find(band)
    FLTMATRIX = sim.FLT.reshape( sim.LCMATRIX.shape )
    bandlist = FLTMATRIX[0,0,0,0,:,0]
    iband = np.where( bandlist==band )[0][0]

    trest = sim.TREST
    ipk = np.argmin( np.abs( trest ) )
    mulist = np.arange(2,28,1)

    outline0 = '# VISIBILITY TIMES FOR %s in %s'%(sim.simname,filter)
    outlinefmt  = '%-3.2f  ' + '%6.2f %6.2f    '*len(mulist)
    outlinefmt1 = '%-3s  ' + '%-12s     '*len(mulist)
    outlinefmt2 = '%-3s  ' + '%-6s %-6s    '*len(mulist)
    outline1 = outlinefmt1 % (
        tuple(np.append( '# ', ['    mu=%i    '%mu for mu in mulist]) ))
    outline2 = outlinefmt2 % (
        tuple(np.append( '# z', [ ['tvis%02i'%mu, 'err%02i'%mu] for mu in mulist]) ))

    fout = open( outfile, 'w')
    print >> fout, outline0
    print >> fout, outline1
    print >> fout, outline2

    if verbose :
        print( outline0 )
        print( outline1 )
        print( outline2 )

    for iz in range(len(sim.z)) :
        tobs = sim.TOBS[ iz ]
        dt = np.mean( np.diff( tobs ) )
        tvismedian = []
        tviserr = []
        for mu in mulist :
            tvislist = []
            for ilp in range(len(sim.LUMIPAR)) :
                mag = sim.LCMATRIX[ ilp, 0, 0, iz, iband, : ] - 2.5*np.log10(mu)
                ivis = np.where( mag <= maglimit  )[0]
                tvis = len( ivis ) * dt
                if mag[-1]<maglimit :
                    # simulated LC did not find the end of the visibility
                    # window, so we extrapolate, requiring a positive slope
                    # (b/c we're working in magnitudes, so a declining LC has
                    # a positive dmag / dt)
                    lcslope = (mag[-1]-mag[-2])/dt
                    if lcslope<=0 : lcslope = 0.05
                    tvis += (maglimit-mag[-1])/lcslope
                tvislist.append( tvis )
                if showplots :
                    pl.plot( tobs, mag, marker=' ', ls='-',
                             label='%.2f'%sim.LUMIPAR[ilp] )
            tvismedian.append( np.median(tvislist) )
            tviserr.append( np.std(tvislist) )
            if showplots :
                pl.legend( loc='upper right' )
                ax = pl.gca()
                ax.axhline( maglimit, ls='--', lw=2, color='0.3')
                ax.set_title( 'z=%.2f mu=%i (return to continue)'%(
                    sim.z[iz], mu ) )
                if not ax.yaxis_inverted() : ax.invert_yaxis()
                pl.draw()
                if verbose :
                    print( 'z=%.2f mu=%i tvis=%.2f +- %.2f'%(
                        sim.z[iz], mu, tvismedian[-1], tviserr[-1] ) )
                userin = raw_input( '(return to continue, pdb to debug)  ')
                if userin == 'pdb' :
                    pdb.set_trace()
                pl.clf()

        outline = outlinefmt % (
            tuple(np.append(z[iz], np.ravel(zip(tvismedian,tviserr)))) )
        if verbose : print( outline )
        print >> fout, outline
    fout.close()

_filter_alpha = { 'F225W':'S','F275W':'T','F336W':'U','F390W':'C',
                'F350LP':'W','F763M':'7', 'F845M':'9',
                'F435W':'B','F475W':'G','F555W':'D','F606W':'V','F625W':'R',
                'F775W':'X','F814W':'I','F850LP':'Z',
                'F125W':'J','F160W':'H','F125W+F160W':'A',
                'F105W':'Y','F110W':'M','F140W':'N',
                'F098M':'L','F127M':'O','F139M':'P','F153M':'Q',
                'G800L':'8','G141':'4','G102':'2','blank':'0',
                'F658N':'6','Blank':'0','BLANK':'0'
                }

def filter2band( filter ):
    """ convert a full ACS or WFC3 filter string into its 
    single-digit alphabetic code.  e.g:
      filter2band('F125W')  ==>  'J'
    """
    import exceptions
    filter = filter.upper()
    if filter in _filter_alpha.keys():
        return( _filter_alpha[filter] )
    else : 
        raise exceptions.RuntimeError(
            "Unknown filter %s"%filter )

def band2filter( alpha ):
    """ convert a single-digit alphabetic bandpass code 
    into a full ACS or WFC3 filter string. e.g:
      band2filter('J') ==> 'F125W'
    """
    import exceptions
    alpha = alpha.upper()
    for filter in _filter_alpha.keys():
        if _filter_alpha[filter] == alpha :
            return( filter )
    else : 
        raise exceptions.RuntimeError(
            "Unknown filter %s"%alpha )


def snapshot_threshold( filter='F140W', etime=30 ):
    """ The Vega magnitude at which we get S/N=10 for a point source,
    in the given filter, for a snapshot visit with total duration of etime
    (in minutes)

     # --------------------------------
     #            etime    S/N=10
     # Filter     [min]  [Vega mag]
       F110W(M)    12      24.7
       F110W(M)    20      25.0
       F110W(M)    30      25.3

       F140W(N)    12      24.0
       F140W(N)    20      24.3
       F140W(N)    30      24.6

    """
    if filter == 'F110W' :
        if etime==12 : return( 24.7 )
        elif etime==20 : return( 25.0 )
        elif etime==30 : return( 25.3 )
    elif filter == 'F140W' :
        if etime==12 : return( 24.0 )
        elif etime==20 : return( 24.3 )
        elif etime==30 : return( 24.6 )


