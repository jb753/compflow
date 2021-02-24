"""
This file contains methods to produce a Turbostream input file given a set of
design parameters. It should not need to be changed and it is not neccesary to
understand what goes on below.
"""
import scipy.optimize, scipy.integrate
from . import compflow
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

# Define a namedtuple to store all information about a stage design
# Note that the reference blade speed is taken at rotor inlet
stage_params = [
    'Yp',  # Stagnation P loss coefficients for stator, rotor
    'Ma',
    'Ma_rel',
    'Al',  # Absolute flow angles
    'Al_rel',  # Relative frame flow angles
    'Lam',  # Reaction
    'Ax_Axin',  # Normalised areas
    'U_sqrt_cpToin',  # Non-dimensional blade speed
]
NonDimStage = namedtuple('NonDimStage', stage_params)


def nondim_stage_from_Al(
                         phi, psi, Al,  # Velocity triangles
                         Ma, ga,  # Compressibility
                         eta,  # Loss
                         Vx_rat=(1., 1.)  # Variation in Vx
                         ):
    """Get N-D geometry, blade speed for set of N-D aerodynamic parameters.

    This compact routine is the only stage design code you will ever need. It
    uses Euler's equation, definitions of the N-D parameters, and standard
    perfect gas compressible flow relations.

    Assumes constant mean radius, angular velocity, hence blade speed.

    From the output of this function, arbitrarily choosing one of Omega or rm
    and providing an inlet stagnation state will completely define the stage in
    dimensional terms.

    TODO - allow for varying the mean radius."""

    # Unpack input
    Ala, Alc = Al

    # First, construct velocity triangles

    # Euler work equation sets tangential velocity upstream of rotor
    # This gives us absolute flow angles everywhere
    tanAlb = np.tan(np.radians(Alc)) * Vx_rat[1] + psi / phi
    Alb = np.degrees(np.arctan(tanAlb))
    Al_all = np.hstack((Ala, Alb, Alc))

    # Get non-dimensional velocities
    Vx_U = np.array([Vx_rat[0], 1., Vx_rat[1]]) * phi
    Vt_U = Vx_U * np.tan(np.radians(Al_all))
    Vtrel_U = Vt_U - 1.
    V_U = np.sqrt(Vx_U**2. + Vt_U**2.)
    Vrel_U = np.sqrt(Vx_U**2. + Vtrel_U**2.)
    Alrel = np.degrees(np.arctan2(Vtrel_U, Vx_U))

    # Use Mach number to get U/cpTo
    V_cpTob = compflow.V_cpTo_from_Ma(Ma, ga)
    U_cpToin = V_cpTob / V_U[1]

    cpToin_Usq = 1./U_cpToin**2
    cpToout_Usq = cpToin_Usq - psi

    # if cpToout_Usq < 0.:
    #     raise ValueError('Negative temperature.')

    # N-D temperatures from input blade Ma and psi
    cpTo_Usq = np.array([cpToin_Usq, cpToin_Usq, cpToout_Usq])

    # Ma nos 
    Ma_all = compflow.Ma_from_V_cpTo(V_U / np.sqrt(cpTo_Usq), ga)
    Ma_rel = Ma_all * Vrel_U / V_U

    # capacity everywhere
    Q = compflow.mcpTo_APo_from_Ma(Ma_all, ga)
    Q_Qin = Q / Q[0]

    # Second, construct annulus line

    # Use polytropic effy to get entropy change
    To_Toin = cpTo_Usq / cpTo_Usq[0]

    # Split evenly between rotor and stator
    Ds_cp = -(1. - 1. / eta) * np.log(To_Toin[-1])
    s_cp = np.hstack((0., 0.5, 1.)) * Ds_cp

    # Convert to stagnation pressures
    Po_Poin = np.exp((ga / (ga - 1.)) * (np.log(To_Toin) + s_cp))

    # Area ratios = span ratios because rm = const
    cosAl = np.cos(np.radians(Al_all))
    Dr_Drin = np.sqrt(To_Toin) / Po_Poin / Q_Qin * cosAl[0] / cosAl

    # Evaluate some other parameters that we might want to target
    T_Toin = To_Toin / compflow.To_T_from_Ma(Ma_all, ga)
    P_Poin = Po_Poin / compflow.Po_P_from_Ma(Ma_all, ga)
    Porel_Poin = P_Poin * compflow.Po_P_from_Ma(Ma_rel, ga)
    Lam = (T_Toin[2] - T_Toin[1]) / (T_Toin[2] - T_Toin[0])

    Yp_vane = (Po_Poin[0] - Po_Poin[1]) / (Po_Poin[0] - P_Poin[1])
    Yp_blade = (Porel_Poin[1] - Porel_Poin[2]) / (Porel_Poin[1] - P_Poin[2])

    # Assemble all of the data into the output object
    dat = NonDimStage(
        Yp=tuple([Yp_vane, Yp_blade]),
        Al=tuple(Al_all),
        Al_rel=tuple(Alrel),
        Ma=tuple(Ma_all),
        Ma_rel=tuple(Ma_rel),
        Ax_Axin=tuple(Dr_Drin),
        Lam=Lam,
        U_sqrt_cpToin=U_cpToin,
    )

    return dat

def nondim_stage_from_Lam(
                          phi, psi, Lam, Alin,  # Velocity triangles
                          Ma, ga,  # Compressibility
                          eta,  # Loss
                          Vx_rat=(1., 1.)
                          ):
    """This function produces an N-D stage design at fixed reaction.

    It is most convenient to evaluate a design in terms of exit flow angle, and
    not reaction. So iterate on the exit flow angle until the target reaction
    is achieved."""

    # Inner function
    def iter_Al(x):
        stg_now = nondim_stage_from_Al(phi, psi, [Alin, x],
                                    Ma, ga, eta, Vx_rat)
        return stg_now.Lam - Lam

    # Solving for Lam in general is tricky
    # Our strategy is to map out a coarse curve first, pick two points
    # that are close to the desired reaction, then secant iterate
    Al_guess = np.linspace(-89.,89.,9)
    with np.errstate(invalid='ignore'):
        Lam_guess = np.array([iter_Al(Ali) for Ali in Al_guess])
    Al_guess = Al_guess[~np.isnan(Lam_guess)]
    Lam_guess = Lam_guess[~np.isnan(Lam_guess)]
    i1, i2 = np.argmax(Lam_guess), np.argmin(Lam_guess)
    Al_guess, Lam_guess = Al_guess[i1:i2], Lam_guess[i1:i2]
    i0 = np.argmin(np.abs(Lam_guess))

    Al_soln = scipy.optimize.newton(iter_Al,x0=Al_guess[i0],x1=Al_guess[i0-1])

    stg_out = nondim_stage_from_Al(phi, psi, [Alin, Al_soln],
                                   Ma, ga, eta, Vx_rat)

    return stg_out

def annulus_line(U_sqrt_cpToin, Ax_Axin, htr, cpToin, Omega):
    """Get dimensional annulus line for given stage design and inlet state.

    Parameters
    ----------
    U_sqrt_cpToin : float
        Non-dimensional blade speed.
    Ax_Axin : array, length 3
        Axial area ratio.
    htr : float
        Hub-to-tip radius ratio.
    cpToin : float
        Inlet specific stagnation enthalpy [J/kg].
    Omega : float
        Shaft angular velocity [rad/s].

    Returns
    -------
    rm : array, length 3
        Mean radius [m].
    Dr : array, length 3
        Annulus span [m].
    """

    # Use non-dimensional blade speed to get U, hence mean radius
    U = U_sqrt_cpToin * np.sqrt(cpToin)
    rm = U / Omega

    # Use hub-to-tip ratio to set span (mdot will therefore float)
    Dr_rm = 2. * (1. - htr) / (1. + htr)
    Dr = rm * Dr_rm * np.array(Ax_Axin) / Ax_Axin[1]

    return rm, Dr

def blade_section(chi):
    """Makes a simple blade geometry from flow angles."""
    # Copy defaults from MEANGEN (Denton)
    tle     = 0.04 # LEADING EDGE THICKNESS/AXIAL CHORD.
    tte     = 0.04 # TRAILING EDGE THICKNESS/AXIAL CHORD.
    tmax   = 0.25 # MAXIMUM THICKNESS/AXIAL CHORD.
    xtmax  = 0.40 # FRACTION OF AXIAL CHORD AT MAXIMUM THICKNESS 
    xmodle   = 0.02 # FRACTION OF AXIAL CHORD OVER WHICH THE LE IS MODIFIED.
    xmodte   = 0.01 # FRACTION OF AXIAL CHORD OVER WHICH THE TE IS MODIFIED.
    tk_typ   = 2.0  # FORM OF BLADE THICKNESS DISTRIBUTION.


    # Camber line
    xhat = 0.5*(1.-np.cos(np.pi * np.linspace(0.,1.,100)))
    tanchi_lim = np.tan(np.radians(chi))
    tanchi = np.interp(xhat,(0.,1.),tanchi_lim)
    yhat = np.insert(scipy.integrate.cumtrapz(tanchi, xhat), 0, 0.)

    power = np.log(0.5)/np.log(xtmax)
    xtrans  = xhat**power

    # Linear thickness component for leading/trailing edges
    tlin = tle + xhat*(tte - tle)

    # Additional component to make up maximum thickness
    tadd = tmax - (tle + xtmax*(tte-tle))

    # Total thickness
    thick  = tlin + tadd*(1.0 - np.abs(xtrans-0.5)**tk_typ/(0.5**tk_typ))

    # Thin the leading and trailing edges
    xhat1 = 1. - xhat
    fac_thin = np.ones_like(xhat)
    fac_thin[xhat<xmodle] = (xhat[xhat<xmodle]/xmodle)**0.3
    fac_thin[xhat1<xmodte] = (xhat1[xhat1<xmodte]/xmodte)**0.3

    # Upper and lower surfaces
    yup = yhat + thick*0.5*fac_thin
    ydown = yhat - thick*0.5*fac_thin

    # # Plot
    # # f,a = plt.subplots()
    # a.plot(xhat,yhat,'k--')
    # a.plot(xhat,yup,'-x')
    # a.plot(xhat,ydown,'-x')
    # a.axis('equal')

    xy = np.vstack((xhat, yup, ydown))

    return xy

def pitch_Zweifel(Z, Al, Ma, Dr, ga, Yp):
    """Calculate pitch-to-chord from Zweifel coefficient."""
    Alr = np.radians(Al)
    Po2_P2 = compflow.Po_P_from_Ma(Ma[1], ga)
    Po2_Po1 = (1. - Yp)/(1.-Yp/Po2_P2)
    Q1 = compflow.mcpTo_APo_from_Ma(Ma[0], ga)
    cosAl = np.cos(Alr)
    V_cpTo_sinAl = compflow.V_cpTo_from_Ma(Ma, ga) * np.sin(Alr)
    s_c = Z * (1.-Po2_Po1/Po2_P2)/Q1/cosAl[0]/(V_cpTo_sinAl[1] - V_cpTo_sinAl[0]) / Dr[0] * np.mean(Dr)
    return np.abs(s_c)

def chord_Re(Re, To1, Po1, Ma, ga, rgas, Yp):
    """Set chord length using Reynolds number."""

    # Get exit stagnation pressure
    Po2_P2 = compflow.Po_P_from_Ma(Ma, ga)
    Po2_Po1 = (1. - Yp)/(1.-Yp/Po2_P2)

    # Exit static state
    P2 = Po1 * Po2_Po1 / Po2_P2
    To2_T2 = compflow.To_T_from_Ma(Ma, ga)
    T2 = To1 / To2_T2
    ro2 = P2 / rgas / T2
    cp = rgas * ga/(ga-1.)
    V2 = compflow.V_cpTo_from_Ma(Ma,ga) * np.sqrt(cp*To1)

    # Get viscosity using 0.62 power approximation
    muref = 1.8e-5
    Tref = 288.0
    mu = muref * (T2/Tref)**0.62

    # Use reynolds number to set chord
    cx = Re * mu / ro2 / V2
    return cx

def meridional_mesh(xc, rm, Dr, c, nr):
    """Generate meridional mesh for a blade row."""

    nxb = 101  # Points along blade chord
    nxb2 = nxb//2  # Points along blade semi-chord

    # Define a cosinusiodal clustering function
    clust = 0.5*(1. - np.cos(np.pi * np.linspace(0.,1.0,nxb)))
    dclust = np.diff(clust)
    dmin, dmax = dclust.min(), dclust.max()

    # Numbers of points in inlet/outlet
    nxu = int((xc[1]-xc[0]-0.5)/dmax)
    nxd = int((xc[3]-xc[2]-0.5)/dmax)

    # Build up the streamwise grid vector
    x = xc[1] + clust  # Blade
    x = np.insert(x[1:], 0, x[0]+(clust[nxb2:]-1.0))  # In front of LE
    x = np.append(x[:-1], x[-1]+(clust[:nxb2]))  # Behind TE
    if nxu>0:
        x = np.insert(x[1:], 0, np.linspace(xc[0],x[0],nxu))  # Inlet
    if nxd>0:
        x = np.append(x[:-1], np.linspace(x[-1],xc[3],nxd))  # Outlet

    # Trim if needed
    x = x[x > xc[0]]
    x = x[x < xc[-1]]

    # Reset endpoints
    x = np.insert(x, 0, xc[0])
    x = np.append(x, xc[-1])

    # Number of radial points
    spf = np.linspace(0.,1.,nr)
    spf2 = np.atleast_2d(spf).T

    # Now annulus lines
    rh = np.interp(x, xc[1:3], rm - Dr/2.)
    rc = np.interp(x, xc[1:3], rm + Dr/2.)
    rh2 = np.atleast_2d(rh)
    rc2 = np.atleast_2d(rc)
    r = spf2*(rc2-rh2) + rh2
    r = r.T

    # Get indices of edges
    i_edge = [np.where(x==xc[i])[0][0] for i in [1,2]]
    i_edge[1] = i_edge[1]+1

    # # Scale by chord
    x = x*c

    # f,a = plt.subplots()
    # a.plot(x,r,'k-')
    # a.axis('equal')

    return x, r, i_edge

def blade_to_blade_mesh(x, r, ii, chi, nrt, s_c):
    """Generate blade section rt, given merid mesh and flow angles."""

    # Define a cosinusiodal clustering function
    clust = 0.5*(1. - np.cos(np.pi * np.linspace(0.,1.0,nrt)))
    clust3 = (clust[...,None,None]).transpose((2,1,0))


    # Chord
    xpass = x[ii[0]:ii[1]] - x[ii[0]]
    c = xpass.ptp()

    nj = np.shape(chi)[1]
    rt = np.empty(np.shape(r) + (nrt,))
    for j in range(nj):

        # Dimensional section
        sect_now = blade_section(chi[:,j])* c
        rt0 = np.interp(xpass,sect_now[0,:],sect_now[1,:])[:,None,None]
        rt1 = np.interp(xpass,sect_now[0,:],sect_now[2,:])[:,None,None]+s_c*c
        rt[ii[0]:ii[1],j,:] = (rt0 + (rt1-rt0)*clust3).squeeze()

    # Now deal with inlet and exit ducts
    # First just extend clustering
    rt[:ii[0],:,:] = rt[ii[0],:,:]
    rt[ii[1]:,:,:] = rt[ii[1]-1,:,:]

    # Set endpoints to a uniform distribution
    unif_rt = np.linspace(0.,1.,nrt)
    unif_rt3 = (unif_rt[...,None,None]).transpose((2,1,0))
    rt[(0,-1),:,:] = rt[(0,-1),:,0][...,None] + (
            rt[(0,-1),:,-1] - rt[(0,-1),:,0])[...,None]*unif_rt3

    # Relax clustering linearly
    lin_x_up = (np.linspace(0.,1.,ii[0]))
    lin_x_up3 = lin_x_up[...,None,None]
    lin_x_dn = np.linspace(1.,0.,len(x)-ii[1])
    lin_x_dn3 = lin_x_dn[...,None,None]
    rt[:ii[0],:,:] = rt[0,:,:][None,...] + (
            rt[ii[0],:,:]-rt[0,:,:])[None,...]*lin_x_up3
    rt[ii[1]:,:,:] = rt[-1,:,:][None,...] + (
            rt[ii[1],:,:]-rt[-1,:,:])[None,...]*lin_x_dn3

    return rt

def make_patch(kind, i, j, k, nxbid=0, nxpid=0, dirs=None)
    # Periodic patches
    p = ts_tstream_type.TstreamPatch()

    p.kind = getattr(ts_tstream_patch_kind, kind)

    p.ist, p.ien = i
    p.jst, p.jen = j
    p.kst, p.ken = k

    p.nxbid = nxbid
    p.nxpid = nxpid

    if dirs is not None:
        p.idir, p.jdir, p.kdir = dirs
    else:
        p.idir, p.jdir, p.kdir = (0,1,2)

    p.nface = 0
    p.nt = 0

    return p

def add_to_grid(g, xin, rin, rtin, ind, bid):
    """From mesh coordinates, add a block with patches to TS grid object"""
    ni, nj, nk = np.shape(rtin)
    rt = rtin+0.
    r = np.repeat(rin[:,:,None],nk,axis=2)
    x = np.tile(xin[:,None,None],(1,nj,nk))
    print(np.shape(x))
    print(np.shape(r))
    print(np.shape(rt))

    # Permute the coordinates into C-style ordering
    # Turbostream is very fussy about this
    xp = np.zeros((nk, nj, ni), np.float32)
    rp = np.zeros((nk, nj, ni), np.float32)
    rtp = np.zeros((nk, nj, ni), np.float32)
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                xp[k,j,i] = x[i,j,k]
                rp[k,j,i] = r[i,j,k]
                rtp[k,j,i] = rt[i,j,k]

    # Generate new block
    b = ts_tstream_type.TstreamBlock()
    b.bid = bid
    b.ni, b.nj, b.nk = ni, nj, nk
    b.np = 0
    b.procid = 0
    b.threadid = 0

    # Add to grid with coordinates
    g.add_block(b)
    for vname, vval in zip(['x','r','rt'],[xp,rp,rtp]):
        g.set_bp(vname, ts_tstream_type.float, bid, vval)

    # Leading and trailing edges
    ile, ite = ii

    # Add periodic patches first

    # Upstream of LE
    periodic_up_1 = make_patch(
        kind='periodic', i=(0,ile+1), j=(0,nj), k=(0,1), dirs=(0,1,6),
        nxbid=bid, nxpid=1,
        )
    periodic_up_2 = make_patch(
        kind='periodic', i=(0,ile+1), j=(0,nj), k=(nk-1,nk), dirs=(0,1,6),
        nxbid=bid, nxpid=0,
        )
    g.add_patch(bid, periodic_up_1)
    g.add_patch(bid, periodic_up_2)

    # Downstream of TE
    periodic_dn_1 = make_patch(
        kind='periodic', i=(ite,ni), j=(0,nj), k=(0,1), dirs=(0,1,6),
        nxbid=bid, nxpid=3,
        )
    periodic_dn_2 = make_patch(
        kind='periodic', i=(ite,ni), j=(0,nj), k=(nk-1,nk), dirs=(0,1,6),
        nxbid=bid, nxpid=2,
        )
    g.add_patch(bid, periodic_dn_1)
    g.add_patch(bid, periodic_dn_2)

    # Slip wall
    slip_j0 = make_patch( kind='slipwall', i=(0,ni), j=(0,1), k=(0,nk) )
    slip_nj = make_patch( kind='slipwall', i=(0,ni), j=(nj-1,nj), k=(0,nk) )
    g.add_patch(bid, slip_j0)
    g.add_patch(bid, slip_nj)

    # Default avs
    for name in ts_tstream_default.av:
        val = ts_tstream_default.av[name]
        if type(val) == type(1):
            g.set_av(name, ts_tstream_type.int, val)
        else:
            g.set_av(name, ts_tstream_type.float, val)

    # Default bvs
    for name in ts_tstream_default.bv:
        for bid in g.get_block_ids():
            val = ts_tstream_default.bv[name]
            if type(val) == type(1):
                g.set_bv(name, ts_tstream_type.int, bid, val)
            else:
                g.set_bv(name, ts_tstream_type.float, bid, val)

def guess_block(g, bid, x, Po, To, Ma, Al, ga, cp, cv, rgas):
    b = g.get_block(bid)
    ni, nj, nk = (b.ni, b.nj, b.nk)
    xb = g.get_bp('x', bid)
    rb = g.get_bp('r', bid)

    # Interpolate guess to block coords
    Pob = np.interp(xb, x, Po)
    Tob = np.interp(xb, x, To)
    Mab = np.interp(xb, x, Ma)
    Alb = np.interp(xb, x, Al)

    # Get velocities
    Vb = cf.from_Ma('V_cpTo', Mab, ga) * np.sqrt(cp * Tob)
    Vxb = Vb * np.cos(np.radians(Alb))
    Vtb = Vb * np.sin(np.radians(Alb))
    Vrb = np.zeros_like(Vb)

    # Static pressure and temperature
    Pb = Pob / cf.from_Ma('Po_P', Mab, ga)
    Tb = Tob / cf.from_Ma('To_T', Mab, ga)

    # Density
    rob = Pb / rgas / Tb

    # Energy
    eb = cv*Tb + (Vb**2.)/2

    # Primary vars
    rovxb = rob * Vxb
    rovrb = rob * Vrb
    rorvtb = rob * rb * Vtb
    roeb = rob * eb

    # Permute the coordinates into C-style ordering
    robp = np.zeros((nk, nj, ni), np.float32)
    rovxbp = np.zeros((nk, nj, ni), np.float32)
    rovrbp = np.zeros((nk, nj, ni), np.float32)
    rorvtbp = np.zeros((nk, nj, ni), np.float32)
    roebp = np.zeros((nk, nj, ni), np.float32)
    for k in range(nk):
        for j in range(nj):
            for i in range(ni):
                robp[k,j,i] = rob[k,j,i]
                rovxbp[k,j,i] = rovxb[k,j,i]
                rovrbp[k,j,i] = rovrb[k,j,i]
                rorvtbp[k,j,i] = rorvtb[k,j,i]
                roebp[k,j,i] = roeb[k,j,i]

    # Apply to grid
    g.set_bp("ro", ts_tstream_type.float, bid, robp)
    g.set_bp("rovx", ts_tstream_type.float, bid, rovxbp)
    g.set_bp("rovr", ts_tstream_type.float, bid, rovrbp)
    g.set_bp("rorvt", ts_tstream_type.float, bid, rorvtbp)
    g.set_bp("roe", ts_tstream_type.float, bid, roebp)


def generate(fname, phi, psi, Lam, Ma, eta, gap_chord, slip_vane, guess_file=None ):
    gamma = 1.4
    rgas = 287.14
    cp = gamma / (gamma - 1.0) * rgas
    cv = cp/gamma

    # Generate a stage in non-dimensional form
    nd_stage = nondim_stage_from_Lam(
            phi=phi,
            psi=psi,
            Lam=Lam,
            Alin=0.0,
            Ma=Ma,
            ga=gamma,
            eta=eta
            )

    # Choose dimensional conditions
    Omega = 50.0 * 2. * np.pi
    htr = 0.99
    Poin = 16e5
    Toin = 1.6e3
    Alin = 0.0
    inlet = State(gamma, rgas, Poin, Toin, Alin)

    Pout = Poin * nd_stage.PR
    # Get radii
    rm, Dr, _, _ = scale_stage(nd_stage, inlet, Omega, htr)

    # Blade sections
    xy_stator = blade_section(nd_stage.chi[:2])
    xy_rotor = blade_section(nd_stage.chi[2:])

    # Get viscosity at the inlet condition using 0.62 power approximation
    muref = 1.8e-5
    Tref = 288.0
    mu = muref * (Toin/Tref)**0.62

    # Use mach to get velocity and density at vane exit
    V_cpTo = cf.from_Ma('V_cpTo',nd_stage.Ma_out[0],gamma)
    Vref = V_cpTo * np.sqrt(cp * Toin)
    rhoo_rho = cf.from_Ma('rhoo_rho',nd_stage.Ma_out[0],gamma)
    roref = Poin / rgas / Toin / rhoo_rho

    # Use reynolds number to set chords
    Re = 3e6
    cx_vane = Re * mu / roref / Vref
    cx_rotor = cx_vane
    cx = np.array((cx_vane, cx_rotor))

    # Put the blades on a common coordinate system
    gap = gap_chord * cx_vane
    xy_stator = xy_stator * cx_vane
    xy_rotor = xy_rotor * cx_rotor
    xy_rotor[0,:] = xy_rotor[0,:] + cx_vane + gap

    # Pitch to chord using Zweifel
    Z = 0.85
    Alr_sta = np.radians(nd_stage.chi[:2])
    Alr_rot = np.radians(nd_stage.chi[2:])
    denom_sta = (np.cos(Alr_sta[1])**2.)*(
            np.tan(Alr_sta[1]) - np.tan(Alr_sta[0]))
    denom_rot = (np.cos(Alr_rot[1])**2.)*(
            np.tan(Alr_rot[1]) - np.tan(Alr_rot[0]))
    s_c = Z / 2. /np.array([denom_sta,-denom_rot])

    # Round to nearest whole number of blades
    nb = np.round(2. * np.pi * rm[0] / (s_c * cx))
    print('nb = ',nb)
    s_c =  2. * np.pi * rm[0] /nb / cx

    xs, rs, rts = row_mesh(xy_stator, rm[:2], Dr[:2], [cx_vane * 3., gap/2.], s_c[0] * cx_vane)

    xr, rr, rtr = row_mesh(xy_rotor, rm[2:], Dr[2:], [gap/2., cx_rotor * 3.], s_c[1] * cx_rotor)

    # Make the mixing planes coincident
    xm = 0.5 * (xs[-1,0,0] + xr[0,0,0])
    xs[-1,:,:] = xm
    xr[0,:,:] = xm

    f, a = plt.subplots()
    a.plot(np.diff(xr[:,0,0]),'-x')
    # Plot out the blade shapes
    f, a = plt.subplots()
    a.plot(xy_stator[0,:],xy_stator[1:,:].T,'-x')
    a.plot(xy_rotor[0,:],xy_rotor[1:,:].T,'-x')
    a.plot(xy_stator[0,:],xy_stator[1:,:].T + s_c[0] * cx_vane,'-x')
    a.plot(xy_rotor[0,:],xy_rotor[1:,:].T + s_c[1] * cx_rotor,'-x')
    a.axis('equal')
    plt.show()

    # Assemble guess
    xg = np.array([-1., xy_stator[0,0],  xy_stator[0,-1],
                        xy_rotor[0,0],  xy_rotor[0,-1], 1.])
    Pog = np.array(nd_stage.Po_Poin) * Poin
    Pog = np.insert(Pog, [0, 1], [Pog[0], Pog[1]])
    Pog = np.append(Pog, Pog[-1])

    Tog = np.array(nd_stage.To_Toin) * Toin
    Tog = np.insert(Tog, [0, 1], [Tog[0], Tog[1]])
    Tog = np.append(Tog, Tog[-1])

    Mag = np.array(nd_stage.Ma_all)
    Mag = np.insert(Mag, [0, 1], [Mag[0], Mag[1]])
    Mag = np.append(Mag, Mag[-1])

    Alg = np.array(nd_stage.Al_all)
    Alg = np.insert(Alg, [0, 1], [Alg[0], Alg[1]])
    Alg = np.append(Alg, Alg[-1])

    print(np.shape(xg))
    print(np.shape(Mag))
    print(Mag)

    # f,a = plt.subplots(4,1)
    # a[0].plot(xg,Pog,'k-x')
    # a[1].plot(xg,Tog,'k-x')
    # a[2].plot(xg,Mag,'k-x')
    # a[3].plot(xg,Alg,'k-x')
    # plt.show()

    g = ts_tstream_grid.TstreamGrid()
    # g = None
    add_to_grid(g, xs, rs, rts, 0, slip_vane)
    add_to_grid(g, xr, rr, rtr, 1)

    for bid in [0,1]:
        guess_block(g, bid, xg, Pog, Tog, Mag, Alg, gamma, cp, cv, rgas)

    # Inlet
    ni, nj, nk = np.shape(xs)
    pin = ts_tstream_type.TstreamPatch()
    pin.kind = ts_tstream_patch_kind.inlet
    pin.bid = 0
    pin.ist = 0
    pin.ien = 1
    pin.jst = 0
    pin.jen = nj
    pin.kst = 0
    pin.ken = nk
    pin.nxbid = 0
    pin.nxpid = 0
    pin.idir = 6
    pin.jdir = 1
    pin.kdir = 2
    pin.pid = g.add_patch(0, pin)

    bid=0
    pid=pin.pid
    print(pin.pid)
    yaw = np.zeros(( nk,  nj), np.float32)
    pitch = np.zeros(( nk,  nj), np.float32)
    pstag = np.zeros(( nk,  nj), np.float32)
    tstag = np.zeros(( nk,  nj), np.float32)

    pstag += Poin
    tstag += Toin
    yaw += 0.0
    pitch += 0.0

    g.set_pp("yaw", ts_tstream_type.float, bid, pid, yaw)
    g.set_pp("pitch", ts_tstream_type.float, bid, pid, pitch)
    g.set_pp("pstag", ts_tstream_type.float, bid, pid, pstag)
    g.set_pp("tstag", ts_tstream_type.float, bid, pid, tstag)
    g.set_pv("rfin", ts_tstream_type.float, bid, pid, 0.5)
    g.set_pv("sfinlet", ts_tstream_type.float, bid, pid, 1.0)

    # Outlet
    ni, nj, nk = np.shape(xr)
    pout = ts_tstream_type.TstreamPatch()
    pout.kind = ts_tstream_patch_kind.outlet
    pout.bid = 1
    pout.ist = ni-1
    pout.ien = ni
    pout.jst = 0
    pout.jen = nj
    pout.kst = 0
    pout.ken = nk
    pout.nxbid = 0
    pout.nxpid = 0
    pout.idir = 6
    pout.jdir = 1
    pout.kdir = 2
    pout.pid = g.add_patch(1, pout)

    g.set_pv("throttle_type", ts_tstream_type.int, pout.bid, pout.pid, 0)
    g.set_pv("ipout", ts_tstream_type.int, pout.bid, pout.pid, 3)
    g.set_pv("pout", ts_tstream_type.float, pout.bid, pout.pid, Pout)

    # mixing upstream
    ni, nj, nk = np.shape(xs)
    pmix1 = ts_tstream_type.TstreamPatch()
    pmix1.kind = ts_tstream_patch_kind.mixing
    pmix1.bid = 0
    pmix1.ist = ni-1
    pmix1.ien = ni
    pmix1.jst = 0
    pmix1.jen = nj
    pmix1.kst = 0
    pmix1.ken = nk
    pmix1.nxbid = 1
    pmix1.nxpid = pout.pid+1
    pmix1.idir = 6
    pmix1.jdir = 1
    pmix1.kdir = 2
    pmix1.pid = g.add_patch(0, pmix1)

    # mixing downstream
    ni, nj, nk = np.shape(xr)
    pmix2 = ts_tstream_type.TstreamPatch()
    pmix2.kind = ts_tstream_patch_kind.mixing
    pmix2.bid = 1
    pmix2.ist = 0
    pmix2.ien = 1
    pmix2.jst = 0
    pmix2.jen = nj
    pmix2.kst = 0
    pmix2.ken = nk
    pmix2.nxbid = 0
    pmix2.nxpid = pmix1.pid
    pmix2.idir = 6
    pmix2.jdir = 1
    pmix2.kdir = 2
    pmix2.pid = g.add_patch(1, pmix2)

    # Mixing length limit
    for bid, nbi in zip([0,1],nb):
        pitchnow = 2. * np.pi * rm[0]  / nbi
        g.set_bv("xllim", ts_tstream_type.float, bid, 0.03*pitchnow)
        g.set_bv("fblade", ts_tstream_type.float, bid, nbi)
        g.set_bv("nblade", ts_tstream_type.float, bid, nbi)

    # other block variables
    for bid in g.get_block_ids():
        g.set_bv("fmgrid", ts_tstream_type.float, bid, 0.2)
        g.set_bv("poisson_fmgrid", ts_tstream_type.float, bid, 0.05)

    g.set_av("restart", ts_tstream_type.int, 1)
    g.set_av("poisson_restart", ts_tstream_type.int, 0)
    g.set_av("poisson_nstep", ts_tstream_type.int, 10000)
    g.set_av("ilos", ts_tstream_type.int, 1)
    g.set_av("nlos", ts_tstream_type.int, 5)
    g.set_av("nstep", ts_tstream_type.int, 16000)
    g.set_av("nchange", ts_tstream_type.int, 1000)
    g.set_av("dampin", ts_tstream_type.float, 5.0)
    g.set_av("sfin", ts_tstream_type.float, 0.5)
    g.set_av("facsecin", ts_tstream_type.float, 0.005)
    g.set_av("cfl", ts_tstream_type.float, 0.4)
    g.set_av("poisson_cfl", ts_tstream_type.float, 0.5)
    g.set_av("fac_stmix", ts_tstream_type.float, 0.0)
    g.set_av("rfmix", ts_tstream_type.float, 0.01)
    g.set_av("viscosity", ts_tstream_type.float, muref)
    g.set_av("viscosity_law", ts_tstream_type.int, 1)

    g.set_av("write_yplus", ts_tstream_type.int, 1)
    # initial guess

    Vguess = Vref * np.cos(Alr_sta[1])
    Pinguess = [Poin, np.mean([Poin,Pout])]
    Poutguess = [np.mean([Poin,Pout]), Pout]
    Toguess = [Toin,Toin*0.8]
    for bid in g.get_block_ids():

        # g.set_bv("ftype", ts_tstream_type.int, bid, 5)
        # g.set_bv("vgridin", ts_tstream_type.float, bid, Vguess)
        # g.set_bv("vgridout", ts_tstream_type.float, bid, Vguess)
        # g.set_bv("pstatin", ts_tstream_type.float, bid, Pinguess[bid])
        # g.set_bv("pstatout", ts_tstream_type.float, bid, Poutguess[bid])
        # g.set_bv("tstagin", ts_tstream_type.float, bid, Toguess[bid])
        # g.set_bv("tstagout", ts_tstream_type.float, bid, Toguess[bid])

        g.set_bv("xllim_free", ts_tstream_type.float, bid, 0.0)
        g.set_bv("free_turb", ts_tstream_type.float, bid, 0.0)

        
    # Rotation
    rpm_rotor = np.array([Omega / 2. / np.pi * 60.]).astype(np.float32)
    for bid, rpm in zip([0,1],[0.,rpm_rotor[0]]):
        print(bid)
        print(rpm)
        g.set_bv("rpm", ts_tstream_type.float, bid, rpm)
        g.set_bv("rpmi1", ts_tstream_type.float, bid, rpm)
        g.set_bv("rpmi2", ts_tstream_type.float, bid, rpm)
        g.set_bv("rpmj1", ts_tstream_type.float, bid, rpm)
        g.set_bv("rpmj2", ts_tstream_type.float, bid, rpm)
        g.set_bv("rpmk1", ts_tstream_type.float, bid, rpm)
        g.set_bv("rpmk2", ts_tstream_type.float, bid, rpm)

    # Load guess from file if provided
    if guess_file is not None:
        print("Loading guess from file...")
        # Load grid
        tsr1 = ts_tstream_reader.TstreamReader()
        g1 = tsr1.read(guess_file)
        # Variables to copy
        varnames = ['ro','rovx','rovr','rorvt','roe']
        # Loop over blocks in new design
        for bid in g.get_block_ids():
            for vi in varnames:
                # Apply flow variables from guess to new grid
                g.set_bp(vi, ts_tstream_type.float, bid, g1.get_bp(vi,bid) )


    g.write_hdf5(fname)

if __name__ == "__main__":

    ga = 1.33
    Toin = 1600.
    Poin = 16.0e5
    rgas = 287.14
    cp = rgas * ga / (ga - 1.)
    htr = 0.6
    Omega = 2.*np.pi*50.
    Z = 0.85
    phi = 0.6
    psi = 1.6
    Lam = 0.5
    Alin=0.
    Ma = 0.75
    eta = 0.95
    Re = 4.0e6

    # Get velocity triangles
    stage = nondim_stage_from_Lam( phi, psi, Lam, Alin, Ma, ga, eta )

    # Get mean radius and spans
    rm, Dr = annulus_line(
        stage.U_sqrt_cpToin, stage.Ax_Axin, htr, cp*Toin, Omega,
        )

    # Calculate pitch-to-chord ratio (const Zweifel)
    s_c = (
        pitch_Zweifel(Z, stage.Al[:2], stage.Ma[:2], Dr[:2], ga, stage.Yp[0]),
        pitch_Zweifel(Z, stage.Al_rel[1:], stage.Ma_rel[1:], Dr[1:], ga, stage.Yp[1])
        )

    # Set chord based on Reynolds number
    c = chord_Re(Re, Toin, Poin, stage.Ma[1], ga, rgas, stage.Yp[0])

    # Number of radial points
    dr_c = 0.05  # Radial grid spacing as fraction of chord
    nr = np.max((int(Dr[1] / (dr_c * c)),4))
    print('nr',nr)

    # Get meridional grids
    vane_xc = [-5.,0.,1., 1.5]
    blade_xc = [1.5,2.,3.,8.]
    vane_x, vane_r, vane_i  = meridional_mesh(vane_xc, rm, Dr[:2], c, nr)
    blade_x, blade_r, blade_i = meridional_mesh(blade_xc, rm, Dr[1:], c, nr)

    # Radial flow angle variations
    vane_Al = np.degrees(np.arctan(
        vane_r[vane_i,:]/rm * np.tan(np.radians(
        np.atleast_2d(stage.Al[:2])).T)))
    blade_Al = np.degrees(np.arctan(
        blade_r[blade_i,:]/rm * np.tan(np.radians(
            np.atleast_2d(stage.Al_rel[1:])).T)))

    nrt = 65
    vane_rt = blade_to_blade_mesh(
        vane_x,
        vane_r,
        vane_i,
        vane_Al,
        nrt,
        s_c[0],
        )
    blade_rt = blade_to_blade_mesh(
        blade_x,
        blade_r,
        blade_i,
        blade_Al,
        nrt,
        s_c[1],
        )

    jp = -1
    f,a = plt.subplots()
    a.plot(vane_x,vane_rt[:,jp,:],'k-')
    a.plot(blade_x,blade_rt[:,jp,:],'r-')
    a.axis('equal')

    add_to_grid(None, vane_x, vane_r, vane_rt, 0, True)
    quit()

    nv = 25
    phi_v = np.linspace(0.4, 1.2, nv)
    psi_v = np.linspace(1.0, 2.4, nv)
    Lam_v = np.linspace(0.4, 0.6, nv)

    phi_g, psi_g, Lam_g = np.meshgrid(phi_v, psi_v, Lam_v)

    Al_soln = []
    invalid = 0
    for phi, psi, Lam in zip(phi_g.flat, psi_g.flat, Lam_g.flat):
        try:
            stage = nondim_stage_from_Lam(
                phi=phi, psi=psi, Lam=Lam, Alin=-10., Ma=1.1, ga=1.33, eta=0.8
                )
            Al_soln.append(stage.Al[-1])
        except RuntimeError:
            invalid += 1

    Al_soln = np.array(Al_soln)
    print(float(invalid)/float(nv**3.))
    print(Al_soln.min(),Al_soln.max())
    plt.hist(Al_soln.flat,100)
    plt.show()
