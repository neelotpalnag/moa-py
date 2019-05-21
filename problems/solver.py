import math
import copy
import random
import numpy as np
import multiprocessing as mp

from optproblems.base import Individual, BoundConstraintsChecker
from optproblems.multiobjective import MultiObjectiveTestProblem
from optproblems.zdt import ZDT1, ZDT2, ZDT3, ZDT4, ZDT6
from optproblems.dtlz import DTLZ1, DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7
from optproblems.wfg import WFG1, WFG2, WFG3, WFG4, WFG5, WFG6, WFG7, WFG8, WFG9


# Solver method
# Args:
# input: nx1 list of float inputs X which translate as individual points in the variable space
# problem : String specifying name of problem
# num_objectives: Integer number of objectives
# **kwargs:
# # multiprocess: True/False
# # maximize: True/Fasle
# # k: int, for WFG Problems

def solver(input, problem, num_objectives, *kwargs):
    mutliprocess = False
    maximize  = False
    k = None

    for key, value in kwargs.items():
        if key == "multiprocess":
            multiprocess = value
            continue
        if key == "maximize":
            maximize = value
            continue
        if key == "k":
            k = value
            continue

    num_variables = len(input)

    # For logging relevant outputs in file
    import logging as l
    l.basicConfig(filename='solver.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s', level=l.INFO)
    l.warning('Starting SOLVER log for %s', problem)

    # Solver should allow only 2 or more than 2 objectives
    if num_objectives<2:
        l.error("SOLVER: Less than two objectives is not allowed. Terminating")
        exit()

    problem = problem.upper()
    Problem = None
    # Nested ifs for problem
    if "ZDT" in problem:
        if problem == "ZDT1":
            Problem = ZDT1(num_variables=num_variables)
        elif problem == "ZDT2":
            Problem = ZDT2(num_variables=num_variables)
        elif problem == "ZDT3":
            Problem = ZDT3(num_variables=num_variables)
        elif problem == "ZDT4":
            Problem = ZDT4(num_variables=num_variables)
        elif problem == "ZDT6":
            Problem = ZDT6(num_variables=num_variables)

    elif "DTLZ" in problem:
        if problem == "DTLZ1":
            Problem = DTLZ1(num_objectives=num_objectives, num_variables=num_variables)
        elif problem == "DTLZ2":
            Problem = DTLZ2(num_objectives=num_objectives, num_variables=num_variables)
        elif problem == "DTLZ3":
            Problem = DTLZ3(num_objectives=num_objectives, num_variables=num_variables)
        elif problem == "DTLZ4":
            Problem = DTLZ4(num_objectives=num_objectives, num_variables=num_variables)
        elif problem == "DTLZ5":
            Problem = DTLZ5(num_objectives=num_objectives, num_variables=num_variables)
        elif problem == "DTLZ6":
            Problem = DTLZ6(num_objectives=num_objectives, num_variables=num_variables)
        elif problem == "DTLZ7":
            Problem = DTLZ7(num_objectives=num_objectives, num_variables=num_variables)

    elif "WFG" in problem:
        # Define the value of k
        # It must hold that k<num_variables
        # Huband et al. recommend k = 4 for two objectives and k = 2 * (m - 1) for m objectives.
        if k == None:
            if num_objectives == 2:
                k = 4
            elif num_objectives>2:
                k = 2*(num_objectives-1)

        if k >= num_variables:
            l.error("WFG: k > Number of Variables is not allowed. Terminating")
            exit()

        if problem == "WFG1":
            Problem = WFG1(num_objectives=num_objectives, num_variables=num_variables, k=k)
        elif problem == "WFG2":
            Problem = WFG2(num_objectives=num_objectives, num_variables=num_variables, k=k)
        elif problem == "WFG3":
            Problem = WFG3(num_objectives=num_objectives, num_variables=num_variables, k=k)
        elif problem == "WFG4":
            Problem = WFG4(num_objectives=num_objectives, num_variables=num_variables, k=k)
        elif problem == "WFG5":
            Problem = WFG5(num_objectives=num_objectives, num_variables=num_variables, k=k)
        elif problem == "WFG6":
            Problem = WFG6(num_objectives=num_objectives, num_variables=num_variables, k=k)
        elif problem == "WFG7":
            Problem = WFG7(num_objectives=num_objectives, num_variables=num_variables, k=k)
        elif problem == "WFG8":
            Problem = WFG8(num_objectives=num_objectives, num_variables=num_variables, k=k)
        elif problem == "WFG9":
            Problem = WFG9(num_objectives=num_objectives, num_variables=num_variables, k=k)

    X = Individual(input)
    Problem.evaluate(X)
    l.info(X.objective_values)
    return(X)
    




#
# z = ZDT1(num_variables=5)
# solution = Individual([0,0.5,1,1,0])
# z.evaluate(solution)
# print(solution.objective_values)
