from bayes_opt import BayesianOptimization
from bayes_opt.event import Events
from bayes_opt.logger import JSONLogger
from bayes_opt.util import load_logs
from functools import partial
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

########
# EX 2 #
########
def run_bayesian_optimization_ex2(ex2, count, low):
    if low:
        file_output_name = 'output_raw_ex2_low_' + str(count) + '.log'
    else:
        file_output_name = 'output_raw_ex2_high_' + str(count) + '.log'
    file = open(file_output_name, "w+")
    file.write('Ex1 Energy; Ex1 Prob; Ex2 Energy; Ex2 Prob; Ex3 Energy; Ex3 Prob; Ex4 Energy; Ex4 Prob; MSE; Chi2; Score\n')
    file.close()

    try:
        os.rmdir('output')
    except:
        pass
    if not os.path.exists('output'):
        os.makedirs('output')

    # Bounds for states
    pbounds = {'ex1': (0., 4.), 'p_ex1': (0., 1.), 'p_ex2': (0., 1.), 'ex3': (0., 4.), 'p_ex3': (0., 1.), 'ex4': (0., 4.), 'p_ex4': (0., 1.)}

    # Partial function
    func_to_optimize = partial(black_box_ex2, ex2=ex2, count=count, low=low)

    # Bayesian Optimizer
    optimizer = BayesianOptimization(f=func_to_optimize, pbounds=pbounds, verbose=2, random_state=42)

    # For saving progress
    if low:
        logger_name = "./logs_ex2_low_" + str(count) + ".json"
    else:
        logger_name = "./logs_ex2_high_" + str(count) + ".json"
    logger = JSONLogger(path=logger_name)
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    # Give it a hint
    optimizer.probe(
        params={"ex1": 0.0, "p_ex1": 0.154, "p_ex2": 0.12, "ex3": 3.6845, "p_ex3": 0.25, "ex4": 3.854, "p_ex4": 0.476},
        lazy=True
    )
    optimizer.maximize(init_points=0, n_iter=0)

    # n_iter: How many steps of bayesian optimization you want to perform
    # init_points: How many steps of random exploration you want to perform
    optimizer.maximize(init_points=100, n_iter=150)

    # print(optimizer.max)

    return

def black_box_ex2(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4, count, low):
    score = run_simulation_ex2(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4, count, low)
    return score

def run_simulation_ex2(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4, count, low):
    global counter
    if low:
        file_output_name = 'output_raw_ex2_low_' + str(count) + '.log'
    else:
        file_output_name = 'output_raw_ex2_high_' + str(count) + '.log'
    file = open(file_output_name, 'a+')

    make_json_file(ex1, p_ex1, ex2, p_ex2, ex3, p_ex3, ex4, p_ex4)

    cmd = '../build/sim config.json -t 32'
    os.system(cmd)

    mse, chi2 = analysis('sim.root', 'output.root', 8)
    score = 100./chi2

    dir_name = 'output/ex2/' + str(counter)

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

def find_bounds_ex2(parameters, base_chi2):
    parameters_ = parameters
    param_value = parameters[2]
    
    param_array = []
    chi2_array = []

    param_array.append(param_value)
    chi2_array.append(base_chi2)

    current_chi2 = base_chi2
    count = 1
    while abs(current_chi2 - base_chi2) < 4:
        current_param = param_value - count*0.01

        run_bayesian_optimization_ex2(current_param, count, True)
        file_ = np.loadtxt('output_raw_ex2_low_' + str(count) + '.log', skiprows=1)
        chi2_ = np.amin(file_[:, 9])

        param_array.append(current_param)
        chi2_array.append(chi2_)

        current_chi2 = chi2_

        if current_chi2 - base_chi2 > 4:
            break

        count = count + 1

    current_chi2 = base_chi2
    count = 1
    while abs(current_chi2 - base_chi2) < 4:
        current_param = param_value + count*0.01

        run_bayesian_optimization_ex2(current_param, count, False)
        file_ = np.loadtxt('output_raw_ex2_high_' + str(count) + '.log', skiprows=1)
        chi2_ = np.amin(file_[:, 9])

        param_array.append(current_param)
        chi2_array.append(chi2_)

        current_chi2 = chi2_

        if current_chi2 - base_chi2 > 4:
            break

        count = count + 1

    bounds_ex2 = open('bounds_ex2.dat', 'w')
    for i in range(len(param_array)):
        bounds_ex2.write("{} {}\n".format(param_array[i], chi2_array[i]))
    bounds_ex2.flush()
    bounds_ex2.close()

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

def find_parameter_errors_fixed():
    num_parameters = 8

    ex1_best = 0.0
    ex2_best = 3.0623
    ex3_best = 3.7088
    ex4_best = 3.8664

    p_ex1_best = 0.1564
    p_ex2_best = 0.1054
    p_ex3_best = 0.2457
    p_ex4_best = 0.5580

    parameters = [ex1_best, p_ex1_best, ex2_best, p_ex2_best, ex3_best, p_ex3_best, ex4_best, p_ex4_best]

    # Get the base chi2
    make_json_file(ex1_best, p_ex1_best, ex2_best, p_ex2_best, ex3_best, p_ex3_best, ex4_best, p_ex4_best)
    cmd = '../build/sim config.json -t 32'
    os.system(cmd)

    base_mse, base_chi2 = analysis('sim.root', 'output.root', num_parameters)

    find_bounds_ex2(parameters, base_chi2) 

if __name__ == '__main__':
    find_parameter_errors_fixed()
