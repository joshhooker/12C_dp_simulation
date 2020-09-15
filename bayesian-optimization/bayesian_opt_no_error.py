from bayes_opt import BayesianOptimization
from bayes_opt.event import Events
from bayes_opt.logger import JSONLogger
from bayes_opt.util import load_logs
import json
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import shutil
import sys

sys.path.insert(0, '../aux_scripts/')

from analysis_no_error import analysis

counter = 0

def run_bayesian_optimization():
    file = open('output_raw.log', "w+")
    file.close()

    # Bounds for states
    pbounds = {'p_ex1': (0., 1.), 'p_ex2': (0., 1.), 'p_ex3': (0., 1.), 'p_ex4': (0., 1.)}

    # Bayesian Optimizer
    # Verbose = 0: Silent
    # Verbose = 1: Prints only when a maximum is observed
    # Verbose = 2
    optimizer = BayesianOptimization(f=black_box, pbounds=pbounds, verbose=2, random_state=42)

    # For saving progress
    logger = JSONLogger(path="./logs.json")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

    # n_iter: How many steps of bayesian optimization you want to perform
    # init_points: How many steps of random exploration you want to perform
    optimizer.maximize(init_points=40, n_iter=40)

    for i, res in enumerate(optimizer.res):
        print("Iteration {}: \n\t{}".format(i, res))

    print(optimizer.max)

def black_box(p_ex1, p_ex2, p_ex3, p_ex4):
    score = run_simulation(p_ex1, p_ex2, p_ex3, p_ex4)
    return score

def run_simulation(p_ex1, p_ex2, p_ex3, p_ex4):
    global counter
    file = open('output_raw.log', 'a+')

    make_json_file(p_ex1, p_ex2, p_ex3, p_ex4)

    cmd = '../build/sim config.json -t ' + str(multiprocessing.cpu_count() - 2)
    os.system(cmd)

    mse = analysis('sim.root', 'output.root')
    score = 100./mse

    dir_name = 'output/' + str(counter)

    try:
        os.mkdir(dir_name)
    except:
        pass

    shutil.move('config.json', str(dir_name) + '/config.json')
    shutil.move('nuclear_states.json', str(dir_name) + '/nuclear_states.json')
    shutil.move('output.root', str(dir_name) + '/output.root')

    try:
        os.remove('sim.root')
    except:
        pass

    file.write('%s %s %s %s %s %s\n' % (p_ex1, p_ex2, p_ex3, p_ex4, mse, score))
    file.close()

    counter += 1

    return score

def make_json_file(p_ex1, p_ex2, p_ex3, p_ex4):
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
    data["beamEnergy"] = 111.22
    data["beamOffset"] = [0, 0]
    data["targetThickness"] = 43.0
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
                {"Energy": 0.000, "Probability": p_ex1},
                {"Energy": 3.089, "Probability": p_ex2},
                {"Energy": 3.685, "Probability": p_ex3},
                {"Energy": 3.854, "Probability": p_ex4}
            ]
        }
    ]

    with open('nuclear_states.json', 'w') as nuclear_file:
        json.dump(nuclear_data, nuclear_file, indent=4)

if __name__ == '__main__':
    run_bayesian_optimization()
