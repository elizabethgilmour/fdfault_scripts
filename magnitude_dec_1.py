import fdfault.analysis
import numpy as np
import matplotlib.pyplot as plt

def load_dat_file(problem, name, datadir):
    slip = fdfault.analysis.output(problem, 'slip', datadir)
    slip.load()
    slip_data = list(slip.fielddata[-1])
    #print(slip_data)
    return slip_data

epsilon = 0.0001
def get_rupture_length(slip_data):
    slip = fdfault.analysis.output(problem, 'slip', datadir)
    
    '''out_data is numpy array
    '''
    lasttimestep = slip_data #get the numpy array of slip for the last time step
    #lasttimestep = lasttimestep.tolist() #convert to a list
    nonzeroind = []
    #print(lasttimestep)
    for i in lasttimestep:
        if i > epsilon:
            nonzeroind.append(lasttimestep.index(i))
    #nonzeroind = np.nonzero(lasttimestep)[0] # gives the indices of all nonzero entries in array form
    if len(nonzeroind) == 0:
        rupture_length = 0.
        return rupture_length
    else:
        #rupture_length = (nonzeroind[-1] - nonzeroind[0])/(0.02) # the distance between the first and last nonzero point divided by number of points per m
        rupture_length = len(nonzeroind)/0.02
        # returns lenght in meters
        #print("nonzero " +  str(nonzeroind))
        #print(len(nonzeroind))
        #print(rupture_length)
        return rupture_length

def calculate_seismic_moment(rupture_length, out_data):
    lasttimestep = out_data #get the numpy array of slip for the last time step
   # print(type(lasttimestep))
    #lasttimestep = lasttimestep.tolist() #convert to a list
    shearmodulus = 3.204e+10 # in pascals
    # gives length of rupture in m
    nonzeroval = []
    for i in lasttimestep:
        if i > epsilon:
            nonzeroval.append(i)
    displacement = np.average(nonzeroval)
    #print(displacement)
    area = rupture_length * rupture_length
    #print(rupture_length)
    
    #displacement = np.average(np.trim_zeros(lasttimestep))
    # average of all nonzero slip in m
    
    seismic_moment = shearmodulus * area * displacement # in Nm
    seismic_moment_dc = (seismic_moment)*10**7 #convert to dyne cm
    moment_magnitude = (0.66 * np.log10(seismic_moment_dc)) - 10.7
    print(moment_magnitude)
    return(moment_magnitude)

datadir = '/Users/mac/Documents/fdfault/data'
problem = 'patch_length'

#out_data = load_dat_file(problem, 'slip', datadir)
#rupture_length = get_rupture_length(out_data)
#moment_magnitude = calculate_seismic_moment(rupture_length, out_data)

HEADERS = ["seed", "momentmagnitude"]


def write_to_csv(data, filename=""):
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


number_of_iterations = 100

name = "patch_length"

problems = [name + str(i).zfill(4) for i in range(90,99)]


for problem in problems:
    try:
        out_data = load_dat_file(problem, "slip", datadir)
        rupture_length = get_rupture_length(out_data)
        moment_magnitude = calculate_seismic_moment(rupture_length, out_data)
        write_to_csv([problem, moment_magnitude], 'patch_lengthmagdec3.txt')
        #problem = str(name + str(i).zfill(4))
        print(problem)
    except ModuleNotFoundError:
        print("file does not exist")
        pass
