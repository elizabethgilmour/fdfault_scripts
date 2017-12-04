import sys
sys.path.append('/Users/pawnoutlet/Documents/fdfault/problems')
import fdfault
import fdfault
import numpy as np
import seistools.rough
import matplotlib.pyplot as plt
import scipy.signal as signal
import seistools.coulomb

# script to create problem and save seed and parameter information to a file

# set up csv file

HEADERS = ["seed", "pore_pressure", "extent", "MUS", "MUD", "sxx", "syy", "sxy", "rms_ratio"]


def write_to_csv(data, filename):
    needs_headers = True
    data = [str(i) for i in data]        
    try:
        with open(filename, "r") as f:
            try:
                if HEADERS[0] in f.readlines()[0]:
                    needs_headers = False
            except IndexError:
                needs_headers = True
    except FileNotFoundError:
        needs_headers = True            
    
    with open(filename, "a") as f:
        if needs_headers:
            f.write(",".join(HEADERS))
            f.write("\n")
            f.write(",".join(data))
        else: 
            f.write(",".join(data))
        f.write("\n")
    return True

# using self similar fault, so hurst = 1.
hurst= np.array([1])


seed = 38

seed_iterations = 20
number_of_iterations = 80
for i in range(number_of_iterations):

	np.random.seed(seed)

	rms_ratio = np.array([np.random.uniform(0.03, 0.001)])
	#var = 0.002

	#p = fdfault.problem('seed' + str(seed_iterations) )
	name = 'patch_length'
	number = str(seed_iterations).zfill(4)
	problem_id = str(name + number)
	p = fdfault.problem(problem_id)

# set rk and fd order

	p.set_rkorder(4)
	p.set_sbporder(4)

	# set time step info

	end_time= 8000
	nt= end_time
	interval_time=200
	p.set_nt(end_time)
	p.set_cfl(0.3)
	p.set_ninfo(interval_time)


	pore_pressure = np.random.uniform(0.01, 1., size = 1)[0]

	# set number of blocks and coordinate information
	nx = 1601
	ny = 401
	p.set_nblocks((1,2,1))
	#p.set_nx_block(([nx], [ny, ny], [1]))
	p.set_nx_block(([nx], [ny, ny], [1]))


	# set block dimensions

	p.set_block_lx((0,0,0),(40.,10.))
	p.set_block_lx((0,1,0),(40.,10.))

	# set block boundary conditions

	p.set_bounds((0,0,0),['absorbing', 'absorbing', 'absorbing', 'none'])
	p.set_bounds((0,1,0),['absorbing', 'absorbing', 'none', 'absorbing'])

	# set block surface
	# need to look at it carefully 

	x = np.linspace(0., 40., nx)
	# _________________________________

	#nx= 3201
	extent = np.random.randint(nx*.3, nx, size = 1)[0]
	extent_km = x[extent]
	#print(extent.min)

	

	#MUS = 0.676
	#MUD = 0.325
	MUD = mud = 0.2+np.random.random()*0.4
	sdrop = (0.6-mud)*np.random.random()+0.2
	MUS = MUD + sdrop
	#sxx=-100.
	#sxy= 42.5
	#syy=-100.
	constant = 0.01
	syy = -np.random.random()*20.-10.
	sxx = syy+(np.random.random()-0.5)*0.5*syy
	sxy = -syy*(mud+constant*sdrop+np.random.random()*sdrop*0.3)

	# sets minimum location (in km) for region in which rupture can occur
	# what should the output be? it's currently an integer but should it be an array?
	min_range= np.where (x>=4)[0][0] 
	# sets maximum location (in km) for region in which rupture can occur
	# I know this should be less than or equal. why isn't it working?
	max_range= np.where (x>= extent_km)[0][0] 
	snst = np.empty(nx)
	normal = np.empty(nx)
	failure_stress = np.empty(nx)
	shear= np.empty(nx)

	#____________________________________


	while True:
		y = 10.*np.ones(nx)+seistools.rough.generate_profile(nx, 40, rms_ratio[0], 10, 1, seed)
		surf = fdfault.curve(nx, 'y', x, y)
		p.set_block_surf((0,0,0), 3, surf)
		p.set_block_surf((0,1,0), 2, surf)
		normx, normy = seistools.rough.generate_normals_2d(x, y, 'y')
		for i in range(nx):
			if i<= extent:
				sxx_p = sxx + pore_pressure
				syy_p = syy + pore_pressure
			else:
				sxx_p = sxx
				syy_p = syy
			# put in for loop for determining whether a value is sxx or sxx + pore_pressure
			sn, st = seistools.coulomb.rotate_xy2nt_2d(sxx_p, sxy, syy_p, (normx[i], normy[i]))
			# normal, shear, snst are lists where each item corresponds to a grid point on the fault
			normal[i]= sn
			shear[i]= st
			snst[i] = abs(st)/abs(sn)
		#	print('stress values', normal[i], shear[i])
		b, a = signal.butter(2, 0.01)
		snflt = signal.filtfilt(b, a, snst)
		maxindex = np.argmax(snflt)        
		if (maxindex >= min_range and maxindex <= max_range):
			pertcenter = surf.get_x(maxindex)
			print('center',pertcenter,'seed', seed, 'index', maxindex)
			break

		seed += 1
	print('normal ' + str(normal))

    
	# Setting the stresses
	
	# set interface type
	p.set_iftype(0,'slipweak')

	# set slip weakening parameters

	p.add_pert(fdfault.swparam('constant', dc = 0.2, mus = MUS, mud = MUD),0)
	p.add_pert(fdfault.swparam('boxcar', x0 = 2. , dx = 2., mus = 10000.),0)
	p.add_pert(fdfault.swparam('boxcar', x0 = 38., dx = 2., mus = 10000.),0)

	# ______________________________________

	# --------------loading stresses through a file
	# add load perturbations
	# add load perturbations through a file

	x_diff=x[2]-x[1]    # finding the difference to know the grid spacing
	patch_length_km = 0.75 # one side patch length. Total will be 2* patch length
	no_grid_points= int(patch_length_km/x_diff)  # this will be the patch length on each side of nucleation center
	friction_coef= MUS
	sn = np.zeros((nx,1))
	s2 = np.zeros((nx,1))
	s3 = np.zeros((nx,1))
	s2[:,0]= shear[:]
	sn[:,0]= normal[:]


	s2[maxindex - no_grid_points : maxindex + no_grid_points,0 ] = abs(normal[maxindex-no_grid_points: maxindex+no_grid_points]) *friction_coef
	s2[maxindex-no_grid_points: maxindex+no_grid_points,0 ] += 0.1 
	print( np.max (s2) )
	p.set_loadfile(0, fdfault.loadfile(nx, 1, sn, s2, s3))

	# add another loadfile for pore pressure

	'''A = np.ones((extent, 1))*pore_pressure
	B = np.zeros((nx - extent, 1))
	sn1 = np.concatenate((A, B), axis = 0)
	s21 = sn1
	start_point = maxindex- no_grid_points # the start grid point of the rupture nucleation zone
	end_point = maxindex + no_grid_points


	s21[start_point : end_point,0 ] = abs(normal[maxindex-no_grid_points: maxindex+no_grid_points]) *friction_coef
	s21[maxindex-no_grid_points: maxindex+no_grid_points,0 ] += 0.1 
	s31 = np.zeros((nx,1))'''
	# adding these to rotated frame of reference
	


	write_to_csv([problem_id, seed, pore_pressure, extent, MUS, MUD, sxx, syy, sxy, rms_ratio[0]], str(name) + ".txt") # this time it creates a new file

	print('pore pressure ' + str(pore_pressure))
	print('pore pressure extent ' +str(extent))
	print('static friction ' + str(MUS))
	print('dynamic friction ' + str(MUD))
	print('sxx ' + str(sxx))
	print('syy ' + str(syy))
	print('shear stress ' + str(sxy))
	print('fault roughness ' + str(rms_ratio))
	print('seed iterations ' + str(seed_iterations))
	print('seed ' + str(seed))
	print('extent in km ' + str(extent_km))
	print(min_range,max_range)

	# -------------------------------------------


	# add output unit


	#p.add_output(fdfault.output('vfault','V',0, end_time-1, interval_time, 0, nx-1, 1, ny, ny, 1, 0, 0, 1))
	#p.add_output(fdfault.output('sfault','S',0, end_time-1, interval_time, 0, nx-1, 1, ny, ny, 1, 0, 0, 1))
	#p.add_output(fdfault.output('snfault','Sn',0, end_time-1, interval_time, 0, nx-1, 1, ny, ny, 1, 0, 0, 1))
	#p.add_output(fdfault.output('sxfault','Sx',0, end_time-1, interval_time, 0, nx-1, 1, ny, ny, 1, 0, 0, 1))
	#p.add_output(fdfault.output('syfault','Sy',0, end_time-1, interval_time, 0, nx-1, 1, ny, ny, 1, 0, 0, 1))
	#p.add_output(fdfault.output('sxxbody','sxx',0, end_time-1, interval_time, 0, 3200, 2, 0, 1601, 2, 0, 0, 1))
	#p.add_output(fdfault.output('sxybody','sxy',0, end_time-1, interval_time, 0, 3200, 2, 0, 1601, 2, 0, 0, 1))
	#p.add_output(fdfault.output('syybody','syy',0, end_time-1, interval_time, 0, 3200, 2, 0, 1601, 2, 0, 0, 1))
	p.add_output(fdfault.output('slip','U',0, end_time-1, interval_time, 0, nx-1, 1, ny, ny, 1, 0, 0, 1))



	
	p.write_input()
	seed= seed+1
	#pore_pressure = pore_pressure + 0.01
	seed_iterations= seed_iterations +1 


