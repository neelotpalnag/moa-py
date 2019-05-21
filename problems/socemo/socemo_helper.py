import matlab.engine

# X ----> NxD array of population
# problem -----> string
# num_objectives -------> int
def socemo_driver(problem, num_objectives, max_evaluations):
    # Use the MATLAB python Engine API to call the SOCEMO function
    engine = matlab.engine.start_matlab()

    # The SOCEMO code accepts the following inputs:
    # problem_file ---> specifying which objective problem to optimize
    # num_objectives
    # max_evaluations ---> max int number of evaluations

    # <TIME START>
    result = engine.socemo(problem, num_objectives, max_evaluations)
    # <TIME END>

    # Terminate matlab engine here
    return(result)
