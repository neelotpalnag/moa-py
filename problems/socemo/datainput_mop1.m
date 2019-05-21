function Data= datainput_mop1(problem, num_objectives)

disp("Solving " + problem + " with " + num_objectives + " objectives:");
Data.dim = 10; %problem dimension variable
Data.xlow=zeros(1,Data.dim); %variable lower bounds
Data.xup=ones(1,Data.dim);  %variable upper bounds
Data.nr_obj = num_objectives; %number of objective functions
Data.objfunction=@(x)my_mofun_p1(x, problem, num_objectives); %handle to black-box objective function evaluation
end %function

function y = my_mofun_p1(x, problem,num_objectives)
%objective function evaluations
global sampledata; %sample data is a global variable and collects all points at whcih we evaluate
%throughout the optimization and the function values; columns 1-Data.dim =
%point, remaining columns = objective function values
y = zeros(1,num_objectives); %initialize objective function value vector (1 row, num of columns = number of objective functions in this example)
%returned objective function values must be in row format

% Generate Objective Values:
ypy = py.solver.solver(x, problem, num_objectives);
ymat = cell(ypy);
y = cellfun(@double, ymat);

sampledata = [sampledata; x(:)',y(1), y(2)];
end
