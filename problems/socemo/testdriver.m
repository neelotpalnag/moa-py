function testdriver
%%testdriver.m does a computationally cheap test run of motv.m. the location of the file testrun.m should
%%be made known to the MATLAB search path
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
close all
datafile = 'datainput_mop4'; %there must be a file with the name "datainput_mop1" which contains all problem information.  
maxeval = 180; %the maximum number of function evaluations we want to allo
setrandseed = 1; %the seed for the random number generator; keep this as the same number if you want to
%replicate a run. If you want to do an ensemble, this number must be changed for each run
socemo(datafile, maxeval, setrandseed) %call the optimization algorithm
load results.mat %load the results (we get the structure Data)
if Data.nr_obj ==2 || Data.nr_obj==3
pareto_plot(Data) %plot the Pareto front (works only for problems with 2 and 3 objectives)
disp(Data.Y_nondom)
end

end%function


% K Values
% DTLZ1 : 5
% DTLZ2 : 10
% DTLZ7 : 20

% MaxEval values for ~380 evals
% ZDT1 : 140
% 


% DTLZ1 : 90
% DTLZ2 : 150
% DTLZ7 : 130


