######################################################################
#
# Simple analysis and plotting script to process problem output
# Made specifically for torus problem. Edited by Leon Chan
# Added packages for 3D volume rendering
#
######################################################################

#import packages#
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os,psutil,sys
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
import multiprocessing as mp

##############################################################################################################################

# Parallelize analysis by spawning several processes using multiprocessing's Pool object
def run_parallel(function,dlist,nthreads):
    pool = mp.Pool(nthreads)
    pool.map_async(function,dlist).get(720000)
    pool.close()
    pool.join()

# Initialize global variables
globalvars_keys = ['PROB','NDIMS','DUMPSDIR','PLOTSDIR']
globalvars = {}
grid = {}

##############################################################################################################################

# Function to generate poloidal (x,z) slice
# Argument must be variable, patch pole (to have x coordinate plotted correctly), averaging in phi option
def xz_slice(var, patch_pole=False, average=False):
	xz_var = np.zeros((2*grid['n1'],grid['n2']))
	if average:
		var = np.mean(var,axis=2)
		for i in range(grid['n1']):
			xz_var[i,:] = var[grid['n1']-1-i,:]
			xz_var[i+grid['n1'],:] = var[i,:]
	else:
		angle = 0.; ind = 0
		for i in range(grid['n1']):
			xz_var[i,:] = var[grid['n1']-1-i,:,ind+grid['n3']//2]
			xz_var[i+grid['n1'],:] = var[i,:,ind]
	if patch_pole:
		xz_var[:,0] = xz_var[:,-1] = 0
	return xz_var

##############################################################################################################################

# Function to generate poloidal (y,z) slice
# Argument must be variable, patch pole (to have y coordinate plotted correctly), averaging in phi option
# Not really called but can include a function call 
def yz_slice(var, patch_pole=False, average=False):
	yz_var = np.zeros((2*grid['n1'],grid['n2']))
	if average:
		var = np.mean(var,axis=2)
		for i in range(grid['n1']):
			yz_var[i,:] = var[grid['n1']-1-i,:]
			yz_var[i+grid['n1'],:] = var[i,:]
	else:
		angle = np.pi/2; ind = np.argmin(abs(grid['phi'][0,0,:]-angle))
		for i in range(grid['n1']):
			yz_var[i,:] = var[grid['n1']-1-i,:,ind+grid['n3']//2]
			yz_var[i+grid['n1'],:] = var[i,:,ind]
	if patch_pole:
		yz_var[:,0] = yz_var[:,-1] = 0
	return yz_var

##############################################################################################################################

# Function to generate toroidal (x,y) slice
# Argument must be variable, averaging in theta option
def xy_slice(var, average=False, patch_phi=False):
    if average:
        xy_var = np.mean(var,axis=1)
    else:
        xy_var = var[:,grid['n2']//2,:]
    #xy_var = np.vstack((xy_var.transpose(),xy_var.transpose()[0])).transpose()
    if patch_phi:
        xy_var[:,0] = xy_var[:,-1] = 0
    return xy_var

##############################################################################################################################

# Function to overlay field lines
# Argument must be axes object, B1, B2 and 'nlines' -> a parameter to account for density of field lines
def plotting_bfield_lines(ax,B1,B2,nlines=20):
    xp = xz_slice(grid['x'], patch_pole=True)
    zp = xz_slice(grid['z'])
    B1_phi_avg = B1.mean(axis=-1) 
    B2_phi_avg = B2.mean(axis=-1)
    AJ_phi = np.zeros([2*grid['n1'],grid['n2']]) 
    for j in range(grid['n2']):
        for i in range(grid['n1']):
            AJ_phi[grid['n1']-1-i,j] = AJ_phi[i+grid['n1'],j] = (np.trapz(grid['gdet'][:i,j,0]*B2_phi_avg[:i,j],dx=grid['dx1']) - np.trapz(grid['gdet'][i,:j,0]*B1_phi_avg[i,:j],dx=grid['dx2']))
    AJ_phi -=AJ_phi.min()
    levels = np.linspace(0,AJ_phi.max(),nlines*2)
    ax.contour(xp, zp, AJ_phi, levels=levels, colors='k')

##############################################################################################################################

# The actual function that computes and plots diagnostics for PROB=torus and NDIMS=3
def analysis_torus3d(dumpval, cmap='plasma', vmin=-5, vmax=0, domain = [-50,50,-50,50], bh=True, shading='gouraud'):

    # initialize #
    plt.clf()

    #print out#
    print("Analyzing {0:04d} dump".format(dumpval))

    #load dump files#
    dfile = h5py.File(os.path.join(globalvars['DUMPSDIR'],'dump_0000{0:04d}.h5'.format(dumpval)),'r')

    #load density, internal energy, velocity, and magnetic field
    rho = dfile['prims'][()][Ellipsis,0]
    uu = np.array(dfile['prims'][()][Ellipsis,1])
    u = np.array(dfile['prims'][()][Ellipsis,2:5])
    B = np.array(dfile['prims'][()][Ellipsis,5:8])

    #gas constant#
    gam = np.array(dfile['header/gam'][()])

    #time#
    t = dfile['t'][()]
    
    #close the file#
    dfile.close()

    #formate the time, 3 sig fig#
    t = "{:.3f}".format(t)

    #take log on density#
    logrho=np.log10(rho)

    #pressure#
    pg = (gam-1)*uu

    #metric tensors#
    gti = grid['gcon'][Ellipsis,0,1:4]
    gij = grid['gcov'][Ellipsis,1:4,1:4]

    #metric related stuff
    beta_i = np.einsum('ijks,ijk->ijks',gti,grid['lapse']**2)

    #lorentz factor#
    qsq = np.einsum('ijky,ijky->ijk',np.einsum('ijkxy,ijkx->ijky',gij,u),u)
    gamma = np.sqrt(1+qsq)

    #velocity#
    ui = u-np.einsum('ijks,ijk->ijks',beta_i,gamma/grid['lapse'])
    ut = gamma/grid['lapse']

    #covariant and contravariant 4-velocity#
    ucon = np.append(ut[Ellipsis,None],ui,axis=3)
    ucov = np.einsum('ijkmn,ijkn->ijkm',grid['gcov'],ucon)

    #magnetic 4-vector#
    bt = np.einsum('ijkm,ijkm->ijk',np.einsum('ijksm,ijks->ijkm',grid['gcov'][Ellipsis,1:4,:],B),ucon)
    bi = (B+np.einsum('ijks,ijk->ijks',ui,bt))/ut[Ellipsis,None]
    bcon = np.append(bt[Ellipsis,None],bi,axis=3)
    bcov = np.einsum('ijkmn,ijkn->ijkm',grid['gcov'],bcon)

    #magnetic energy
    bsq = np.einsum('ijkm,ijkm->ijk',bcon,bcov)

    #magnetization
    betainv = 0.5*bsq/pg
    logbetainv = np.log10(betainv) #Log($\\beta^{-1}#

    #x, z slice#
    xp = xz_slice(grid['x'], patch_pole=True)
    zp = xz_slice(grid['z'])
    rhop = xz_slice(logrho)
    betainvp = xz_slice(logbetainv)

    #x, y slice#
    xt = xy_slice(grid['x'])
    yt = xy_slice(grid['y'],patch_phi=True)
    rhot = xy_slice(logrho)
    betainvt = xy_slice(logbetainv)

    #plot density, x-z#
    fig, ax = plt.subplots()
    rhopolplot = ax.pcolormesh(xp, zp, rhop, cmap=cmap, vmin=vmin, vmax=vmax, shading=shading)
    plotting_bfield_lines(ax,B[Ellipsis,0],B[Ellipsis,1],nlines=10)
    ax.set_xlabel('$x (GM/c^2)$',size=15)
    ax.set_ylabel('$z (GM/c^2)$',size=15)
    ax.set_xlim(domain[:2])
    ax.set_ylim(domain[2:])
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_title('Log($\\rho$)',fontsize='large',size=15)
    if bh:
            circle = plt.Circle((0,0),grid['rEH'],color='k')
            ax.add_artist(circle)
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(rhopolplot, cax=cax)
    cbar.ax.tick_params(labelsize=15)
    plt.tight_layout()
    plt.savefig(os.path.join(globalvars['PLOTSDIR'],'{}_logrhoxz_plot_{:04d}.png'.format(globalvars['PROB'],dumpval)))
    plt.clf()
    plt.close()

    #plot density, x-y#
    fig, ax = plt.subplots()
    rhopolplot = ax.pcolormesh(xt, yt, rhot, cmap=cmap, vmin=vmin, vmax=vmax, shading=shading)
    ax.set_xlabel('$x (GM/c^2)$',size=15)
    ax.set_ylabel('$y (GM/c^2)$',size=15)
    ax.set_xlim(domain[:2])
    ax.set_ylim(domain[2:])
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_title('Log($\\rho$)',fontsize='large',size=15)
    if bh:
            circle = plt.Circle((0,0),grid['rEH'],color='k')
            ax.add_artist(circle)
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(rhopolplot, cax=cax)
    cbar.ax.tick_params(labelsize=15)
    plt.tight_layout()
    plt.savefig(os.path.join(globalvars['PLOTSDIR'],'{}_logrhoxy_plot_{:04d}.png'.format(globalvars['PROB'],dumpval)))
    plt.clf()
    plt.close()

##############################################################################################################################

# main(): Reads param file, writes grid dict and calls analysis function
if __name__=="__main__":

    # my patch, the dump file is fixed #
    globalvars['PROB'] = 'torus'
    globalvars['NDIMS'] = '3'
    globalvars['DUMPSDIR'] = './dumps'
    globalvars['PLOTSDIR'] = './figure'

    # Calculating total dump files
    dstart = int(sorted(os.listdir(globalvars['DUMPSDIR']))[0][-7:-3])
    dend = int(sorted(list(filter(lambda dump: 'dump' in dump,os.listdir(globalvars['DUMPSDIR']))))[-1][-7:-3])
    dlist = range(dstart,dend+1)
    Ndumps = dend-dstart+1

    # Setting grid dict
    gfile = h5py.File(os.path.join(globalvars['DUMPSDIR'],'grid.h5'),'r')
    dfile = h5py.File(os.path.join(globalvars['DUMPSDIR'],'dump_0000{0:04d}.h5'.format(dstart)),'r')

    #load number of grid points#
    grid['n1'] = dfile['/header/n1'][()]; grid['n2'] = dfile['/header/n2'][()]; grid['n3'] = dfile['/header/n3'][()]

    #load dx#
    grid['dx1'] = dfile['/header/geom/dx1'][()]; grid['dx2'] = dfile['/header/geom/dx2'][()]; grid['dx3'] = dfile['/header/geom/dx3'][()]

    #load grid starting point#
    grid['startx1'] = dfile['header/geom/startx1'][()]; grid['startx2'] = dfile['header/geom/startx2'][()]; grid['startx3'] = dfile['header/geom/startx3'][()]

    #load metric tensors#
    grid['metric'] = dfile['header/metric'][()].decode('UTF-8')

    #load metric variables, for kerr schild
    if grid['metric']=='MKS' or grid['metric']=='MMKS':
        try:
            grid['a'] = dfile['header/geom/mks/a'][()]
        except KeyError:
            grid['a'] = dfile['header/geom/mmks/a'][()]
        try:
            grid['rEH'] = dfile['header/geom/mks/Reh'][()]
        except KeyError:
            pass
        try:
            grid['rEH'] = dfile['header/geom/mks/r_eh'][()]
        except KeyError:
            pass
        try:
            grid['rEH'] = dfile['header/geom/mmks/Reh'][()]
        except KeyError:
            pass
        try:
            grid['rEH'] = dfile['header/geom/mmks/r_eh'][()]
        except KeyError:
            pass
        try:
            grid['hslope'] = dfile['header/geom/mks/hslope'][()]
        except KeyError:
            grid['hslope'] = dfile['header/geom/mmks/hslope'][()]
    
    #load metric variables, for modified kerr schild
    if grid['metric']=='MMKS':
        grid['mks_smooth'] = dfile['header/geom/mmks/mks_smooth'][()]
        grid['poly_alpha'] = dfile['header/geom/mmks/poly_alpha'][()]
        grid['poly_xt'] = dfile['header/geom/mmks/poly_xt'][()]
        grid['D'] = (np.pi*grid['poly_xt']**grid['poly_alpha'])/(2*grid['poly_xt']**grid['poly_alpha']+(2/(1+grid['poly_alpha'])))

    #load grid parameters from grid files, x1, x2, x3 #
    grid['x1'] = gfile['X1'][()]; grid['x2'] = gfile['X2'][()]; grid['x3'] = gfile['X3'][()]

    #r, theta, phi#
    grid['r'] = gfile['r'][()]; grid['th'] = gfile['th'][()]; grid['phi'] = gfile['phi'][()]

    #x, y, z computed from r, theta, phi#
    grid['x'] = gfile['X'][()]; grid['y'] = gfile['Y'][()]; grid['z'] = gfile['Z'][()]

    # convariant and contravariant tensors#
    grid['gcov'] = gfile['gcov'][()]; grid['gcon'] = gfile['gcon'][()]

    #determinat#
    grid['gdet'] = gfile['gdet'][()]

    #lapse functions#
    grid['lapse'] = gfile['lapse'][()]

    #close all of them#
    dfile.close()
    gfile.close()

    #parallel processing stuff
    ncores = psutil.cpu_count(logical=True)
    pad = 0.25
    nthreads = int(ncores*pad); print("Number of threads: {0:03d}".format(nthreads))

    # Calling analysis function for torus3d
    if globalvars['PROB']=='torus':
        run_parallel(analysis_torus3d,dlist,nthreads)
    else:
        exit()