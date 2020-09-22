from bayes_opt import BayesianOptimization
from bayes_opt.event import Events
from bayes_opt.logger import JSONLogger
from bayes_opt.util import load_logs
import json
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
from scipy.interpolate import CubicSpline
import shutil
import sys

sys.path.insert(0, '../aux_scripts/')

from analysis_no_error import analysis

counter = 0

def run_bayesian_optimization():
    file = open('output_raw.log', "w+")
    file.write('Ex1 Energy; Ex1 Prob; Ex2 Energy; Ex2 Prob; Ex3 Energy; Ex3 Prob; Ex4 Energy; Ex4 Prob; MSE; Chi2; Score\n')
    file.close()

    try:
        os.rmdir('output')
    except:
        pass
    if not os.path.exists('output'):
        os.makedirs('output')

    # Bounds for states
    pbounds = {'ex1': (0., 4.), 'p_ex1': (0., 1.), 'ex2': (0., 4.), 'p_ex2': (0., 1.), 'ex3': (0., 4.), 'p_ex3': (0., 1.), 'ex4': (0., 4.), 'p_ex4': (0., 1.)}

    # Bayesian Optimizer
    # Verbose = 0: Silent
    # Verbose = 1: Prints only when a maximum is observed
    # Verbose = 2
    optimizer = BayesianOptimization(f=black_box, pbounds=pbounds, verbose=2, random_state=42)

    # For saving progress
    logger = JSONLogger(path="./logs.json")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    # Give it a hint
    optimizer.probe(
        params=[0.0, 0.154, 3.089, 0.12, 3.6845, 0.25, 3.854, 0.476],
        lazy=True
    )
<<<<<<< HEAD

    # n_iter: How many steps of bayesian optimization you want to perform
    # init_points: How many steps of random exploration you want to perform
    optimizer.maximize(init_points=150, n_iter=150)
=======
    optimizer.maximize(init_points=0, n_iter=0)

    # n_iter: How many steps of bayesian optimization you want to perform
    # init_points: How many steps of random exploration you want to perform
    optimizer.maximize(init_points=100, n_iter=150)
>>>>>>> Added bayesopt hint

    print(optimizer.max)

    return optimizer.max, len(pbounds)

def black_box(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4):
    score = run_simulation(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4)
    return score

def run_simulation(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4):
    global counter
    file = open('output_raw.log', 'a+')

    make_json_file(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4)

    cmd = '../build/sim config.json -t 10'
    os.system(cmd)

    mse, chi2 = analysis('sim.root', 'output.root', 8)
    score = 100./mse

    dir_name = 'output/' + str(counter)

    try:
        os.mkdir(dir_name)
    except:
        pass

    shutil.move('config.json', str(dir_name) + '/config.json')
    shutil.move('nuclear_states.json', str(dir_name) + '/nuclear_states.json')
    shutil.move('sim.root', str(dir_name) + '/sim.root')

    try:
        os.remove('sim.root')
    except:
        pass

    print(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4, mse, chi2, score)
    file.write('%s %s %s %s %s %s %s %s %s %s %s\n' % (ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4, mse, chi2, score))
    file.close()

    counter += 1

    return score

def make_json_file(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4):
    try:
        os.remove('config.json')
        os.remove('nuclear_states.json')
    except:
        pass

    data = {}
    data["interactive"] = False
    data["macroName"] = "run.mac"
    data["outputFile"] = "sim"
    data["beamIon[Z,A]"] = [6, 12]
    data["targetIon[Z,A]"] = [1, 2]
    data["ejectileIon[Z,A]"] = [1, 1]
    data["beamEnergy"] = 111.4
    data["beamOffset"] = [0, 0]
    data["targetThickness"] = 52.0
    data["yuEnergyFWHM"] = 0.035
    data["yuThreshold"] = 0.2
    data["use13BeStates"] = False
    data["stateEnergyFWHM"] = 0.1
    data["backgroundAmount"] = 0.0
    data["energyQvalue"] = False

    with open('config.json', 'w') as outfile:
        json.dump(data, outfile, indent=4)

    nuclear_data = {}
    nuclear_data["Isotopes"] = [
        {
            "Name": "13C",
            "ZA": [6, 13],
            "States": [
                {"Energy": ex1, "Probability": p_ex1},
                {"Energy": ex2, "Probability": p_ex2},
                {"Energy": ex3, "Probability": p_ex3},
                {"Energy": ex4, "Probability": p_ex4}
            ]
        }
    ]

    with open('nuclear_states.json', 'w') as nuclear_file:
        json.dump(nuclear_data, nuclear_file, indent=4)

def find_bounds(parameters, number, base_chi2):
    parameters_ = parameters
    param_value = parameters[number]
    print(base_chi2, len(parameters))

    param_array = []
    chi2_array = []

    param_array.append(param_value)
    chi2_array.append(base_chi2)

    current_chi2 = base_chi2
    count = 1
    while current_chi2 - base_chi2 < 2:
        current_param = param_value - count*0.01
        if current_param < 0:
            break
        parameters_[number] = current_param
        make_json_file(parameters_[0], parameters_[1], parameters_[2], parameters_[3],
                       parameters_[4], parameters_[5], parameters_[6], parameters_[7])
        cmd = '../build/sim config.json -t 10'
        os.system(cmd)
        mse, chi2 = analysis('sim.root', 'output.root', len(parameters))
        param_array.append(current_param)
        chi2_array.append(chi2)
        print(current_param, chi2)
        current_chi2 = chi2
        if current_chi2 - base_chi2 > 2:
            break
        count = count + 1

    current_chi2 = base_chi2
    count = 1
    while current_chi2 - base_chi2 < 2:
        current_param = param_value + count*0.01
        parameters_[number] = current_param
        make_json_file(parameters_[0], parameters_[1], parameters_[2], parameters_[3],
                       parameters_[4], parameters_[5], parameters_[6], parameters_[7])
        cmd = '../build/sim config.json -t 10'
        os.system(cmd)
        mse, chi2 = analysis('sim.root', 'output.root', len(parameters))
        param_array.append(current_param)
        chi2_array.append(chi2)
        print(current_param, chi2)
        current_chi2 = chi2
        if current_chi2 - base_chi2 > 2:
            break
        count = count + 1

    param_array, chi2_array = zip(*sorted(zip(param_array, chi2_array)))

    return param_array, chi2_array


def find_parameter_errors(max_parameters, num_parameters):
    ex1_best = max_parameters['params']['ex1']
    ex2_best = max_parameters['params']['ex2']
    ex3_best = max_parameters['params']['ex3']
    ex4_best = max_parameters['params']['ex4']

    p_ex1_best = max_parameters['params']['p_ex1']
    p_ex2_best = max_parameters['params']['p_ex2']
    p_ex3_best = max_parameters['params']['p_ex3']
    p_ex4_best = max_parameters['params']['p_ex4']

    parameters = [ex1_best, p_ex1_best, ex2_best, p_ex2_best, ex3_best, p_ex3_best, ex4_best, p_ex4_best]

    # Get the base chi2
    make_json_file(ex1_best, p_ex1_best, ex2_best, p_ex2_best, ex3_best, p_ex3_best, ex4_best, p_ex4_best)
    cmd = '../build/sim config.json -t 10'
    os.system(cmd)

    base_mse, base_chi2 = analysis('sim.root', 'output.root', num_parameters)

    ex1_energy, ex1_chi2 = find_bounds(parameters, 0, base_chi2)
    ex2_energy, ex2_chi2 = find_bounds(parameters, 2, base_chi2)
    ex3_energy, ex3_chi2 = find_bounds(parameters, 4, base_chi2)
    ex4_energy, ex4_chi2 = find_bounds(parameters, 6, base_chi2)

    p_ex1_energy, p_ex1_chi2 = find_bounds(parameters, 1, base_chi2)
    p_ex2_energy, p_ex2_chi2 = find_bounds(parameters, 3, base_chi2)
    p_ex3_energy, p_ex3_chi2 = find_bounds(parameters, 5, base_chi2)
    p_ex4_energy, p_ex4_chi2 = find_bounds(parameters, 7, base_chi2)

    file_ex1 = open('ex1_parameter_estimation.dat', "w+")
    for i in range(len(ex1_energy)):
        file_ex1.write('%s %s\n' % (ex1_energy[i], ex1_chi2[i]))
    file_ex1.close()

    file_ex2 = open('ex2_parameter_estimation.dat', "w+")
    for i in range(len(ex2_energy)):
        file_ex2.write('%s %s\n' % (ex2_energy[i], ex2_chi2[i]))
    file_ex2.close()

    file_ex3 = open('ex3_parameter_estimation.dat', "w+")
    for i in range(len(ex3_energy)):
        file_ex3.write('%s %s\n' % (ex3_energy[i], ex3_chi2[i]))
    file_ex3.close()

    file_ex4 = open('ex4_parameter_estimation.dat', "w+")
    for i in range(len(ex4_energy)):
        file_ex4.write('%s %s\n' % (ex4_energy[i], ex4_chi2[i]))
    file_ex4.close()

    file_p_ex1 = open('p_ex1_parameter_estimation.dat', "w+")
    for i in range(len(p_ex1_energy)):
        file_p_ex1.write('%s %s\n' % (p_ex1_energy[i], p_ex1_chi2[i]))
    file_p_ex1.close()

    file_p_ex2 = open('p_ex2_parameter_estimation.dat', "w+")
    for i in range(len(p_ex2_energy)):
        file_p_ex2.write('%s %s\n' % (p_ex2_energy[i], p_ex2_chi2[i]))
    file_p_ex2.close()

    file_p_ex3 = open('p_ex3_parameter_estimation.dat', "w+")
    for i in range(len(p_ex3_energy)):
        file_p_ex3.write('%s %s\n' % (p_ex3_energy[i], p_ex3_chi2[i]))
    file_p_ex3.close()

    file_p_ex4 = open('p_ex4_parameter_estimation.dat', "w+")
    for i in range(len(p_ex4_energy)):
        file_p_ex4.write('%s %s\n' % (p_ex4_energy[i], p_ex4_chi2[i]))
    file_p_ex4.close()

def find_parameter_errors_fixed():
    num_parameters = 8

    ex1_best = 0.0
    ex2_best = 3.565
    ex3_best = 3.885
    ex4_best = 3.825

    p_ex1_best = 0.2699
    p_ex2_best = 0.6542
    p_ex3_best = 0.9627
    p_ex4_best = 0.0622

    parameters = [ex1_best, p_ex1_best, ex2_best, p_ex2_best, ex3_best, p_ex3_best, ex4_best, p_ex4_best]

    # Get the base chi2
    make_json_file(ex1_best, p_ex1_best, ex2_best, p_ex2_best, ex3_best, p_ex3_best, ex4_best, p_ex4_best)
    cmd = '../build/sim config.json -t 10'
    os.system(cmd)

    base_mse, base_chi2 = analysis('sim.root', 'output.root', num_parameters)

    ex1_energy, ex1_chi2 = find_bounds(parameters, 0, base_chi2)
    ex2_energy, ex2_chi2 = find_bounds(parameters, 2, base_chi2)
    ex3_energy, ex3_chi2 = find_bounds(parameters, 4, base_chi2)
    ex4_energy, ex4_chi2 = find_bounds(parameters, 6, base_chi2)

    p_ex1_energy, p_ex1_chi2 = find_bounds(parameters, 1, base_chi2)
    p_ex2_energy, p_ex2_chi2 = find_bounds(parameters, 3, base_chi2)
    p_ex3_energy, p_ex3_chi2 = find_bounds(parameters, 5, base_chi2)
    p_ex4_energy, p_ex4_chi2 = find_bounds(parameters, 7, base_chi2)

    file_ex1 = open('ex1_parameter_estimation.dat', "w+")
    for i in range(len(ex1_energy)):
        file_ex1.write('%s %s\n' % (ex1_energy[i], ex1_chi2[i]))
    file_ex1.close()

    file_ex2 = open('ex2_parameter_estimation.dat', "w+")
    for i in range(len(ex2_energy)):
        file_ex2.write('%s %s\n' % (ex2_energy[i], ex2_chi2[i]))
    file_ex2.close()

    file_ex3 = open('ex3_parameter_estimation.dat', "w+")
    for i in range(len(ex3_energy)):
        file_ex3.write('%s %s\n' % (ex3_energy[i], ex3_chi2[i]))
    file_ex3.close()

    file_ex4 = open('ex4_parameter_estimation.dat', "w+")
    for i in range(len(ex4_energy)):
        file_ex4.write('%s %s\n' % (ex4_energy[i], ex4_chi2[i]))
    file_ex4.close()

    file_p_ex1 = open('p_ex1_parameter_estimation.dat', "w+")
    for i in range(len(p_ex1_energy)):
        file_p_ex1.write('%s %s\n' % (p_ex1_energy[i], p_ex1_chi2[i]))
    file_p_ex1.close()

    file_p_ex2 = open('p_ex2_parameter_estimation.dat', "w+")
    for i in range(len(p_ex2_energy)):
        file_p_ex2.write('%s %s\n' % (p_ex2_energy[i], p_ex2_chi2[i]))
    file_p_ex2.close()

    file_p_ex3 = open('p_ex3_parameter_estimation.dat', "w+")
    for i in range(len(p_ex3_energy)):
        file_p_ex3.write('%s %s\n' % (p_ex3_energy[i], p_ex3_chi2[i]))
    file_p_ex3.close()

    file_p_ex4 = open('p_ex4_parameter_estimation.dat', "w+")
    for i in range(len(p_ex4_energy)):
        file_p_ex4.write('%s %s\n' % (p_ex4_energy[i], p_ex4_chi2[i]))
    file_p_ex4.close()

if __name__ == '__main__':
    max_parameters, num_parameters = run_bayesian_optimization()
    find_parameter_errors(max_parameters, num_parameters)

    # find_parameter_errors_fixed()
