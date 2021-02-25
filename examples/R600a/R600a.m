%% RankineLab - A MATLAB program to optimize Rankine cycles
% Authors: Roberto Agromayor
% Date: Spring 2021


%% Initialize the program
% Clear all variables and close all figures
clear all
close all
clc

% Add the path to the 'source_code' directory
% The function genpath() is very convenient
addpath(genpath('../../source_code'))

% Set a project name (mfilename corresponds to the script name)
project_name = mfilename;

% Create a directory to store the figures and results
results_path = fullfile(pwd, [project_name, '_results']);
if exist(results_path, 'dir') ~= 7    % 7 is a MATLAB's convention
    mkdir(results_path)
else
    rmdir(results_path, 's')
    mkdir(results_path)
end


%% Define the thermodynamic specifications and cycle parameters
% Ambient state
p_0 = 101325;                                                              % Ambient pressure
T_0 = 10+273.15;                                                           % Ambient temperature

% Define the working fluids
heating_fluid = 'HEOS::Air';
working_fluid = 'HEOS::Isobutane';
cooling_fluid = 'HEOS::Water';

% Some guidance to select the fluid name and EoS
%{

You can set the back-end used to compute the thermodynamic properties by
precedding the fluid name by 'HEOS::' (for the Coolprop back-end) or by
'REFPROP::' (for the REFPROP backend).

You can check the list of the CoolProp back-end fluid names here:
   http://www.coolprop.org/fluid_properties/PurePseudoPure.html

Note that the fluid names for REFPROP may be different

The equations of state may fail for some thermodynamic states and break
down the optimization progress (oh, the horror...)

This is depedent on the fluid and the EoS back-end (REFPROP/CoolProp)

My recommendation is to try the CoolProp back-end first
   1) CoolProp EoS are (much) faster, but can fail sometimes
   2) REFPROP EoS are slower, but they are also a bit more robust

The worst kind of error that you can get while executing the code is
that the equations of state fail. 

There is not much to do about it when it happens because it depends on the
CoolProp/REFPROP libraries.

If you encounter this problem you can try to start the optimization from a
different initial guess, change the bounds/constraints of the optimization
or, as a last measure, change the EoS back-to REFPROP.

If all of these fail, you are out of luck...

%}

% Define the thermodynamic specifications of the heat source
m_h = 1.00;                                                                % Heat source mass flow rate (kg/s). This variable "scales-up" the problem
T_1 = 150+273.15;                                                          % Heat source inlet temperature (K)
p_1 = p_0;                                                                 % Heat source inlet pressure (Pa)
p_3 = p_0;                                                                 % Heat source outlet pressure (Pa)
T_3min = 80+273.15;                                                        % Minimum temperature at the outlet of the heat source (K)
T_3max = T_1;                                                              % Maximum temperature at the outlet of the heat source (K)

% Define the thermodynamic specifications of the heat sink
T_10 = T_0;                                                                % Heat sink inlet temperature (K)
p_10 = p_0;                                                                % Heat sink inlet pressure (Pa)
p_12 = p_0;                                                                % Heat sink outlet pressure (Pa)
T_12min = T_0+10;                                                          % Minimum temperature at the outlet of the heat sink (K)
T_12max = T_0+20;                                                          % Maximum temperature at the outlet of the heat sink (K)

% Define pressure drops in the heat exchangers (% of the inlet pressure)
dp_h_evap = 0.02;                                                          % Pressure drop in the heating fluid side of the evaporator (relative with respect to inlet)
dp_c_evap = 0.02;                                                          % Pressure drop in the working fluid side of the evaporator (relative with respect to inlet)
dp_h_cond = 0.02;                                                          % Pressure drop in the working fluid side of the condenser (relative with respect to inlet)
dp_c_cond = 0.02;                                                          % Pressure drop in the cooling fluid side of the condenser (relative with respect to inlet)
dp_h_rec  = 0.00;                                                          % Pressure drop in the hot side of the recuperator (relative with respect to inlet)
dp_c_rec  = 0.00;                                                          % Pressure drop in the cold side of the recuperator (relative with respect to inlet)

% Define the number of discretizations of the heat exchangers
% More discretizations: more accurate pinch point and more execution time
N_evap = 50;                                                               % Choose an integer
N_cond = 50;                                                               % Choose an integer
N_rec  = 50;                                                               % Choose an integer

% Define the efficiency of the pumps
eta_pump_h = 0.80;                                                         % Isentropic efficiency of the heat source pump
eta_pump_f = 0.80;                                                         % Isentropic efficiency of the working fluid pump
eta_pump_c = 0.80;                                                         % Isentropic efficiency of the heat sink pump

% Define the efficiency of the expander
eta_expander = 0.85;                                                       % Efficiency of the expander

% Turbomachinery efficiency definition
pump_efficiency_definition = 'isentropic';                                 % Definition of the pump efficiency: 'isentropic' or 'polytropic'
expander_efficiency_definition = 'isentropic';                             % Definition of the expander efficiency: 'isentropic' or 'polytropic'

% More details about the efficiency definition
%{

If you set efficiency_definition='isentropic' the efficiency of the
working fluid pump and expander are defined as:
    
   eta_pump = (h_out_s-h_in)/(h_out-h_in)
   eta_expander = (h_in-h_out)/(h_in-h_out_s)

and the turbomachinery computations only involve the inlet and exit states

If you set efficiency_definition='polytropic' the efficiency of the working
fluid pump and expander are defined as:

   eta_pump = dh_s/dh = 1/d*dp/dh
   eta_expander = dh/dh_s = d*dh/dp

and the turbomachinery computations are performed by solving the previous
ordinary differential equations (ODE) along the compression or expansion

%}


%% Define extra fluid parameters
% Used to scale degrees of freedom (no need to do changes in most cases)
% Define the critical properties
p_crit = prop_calculation('P_CRITICAL',working_fluid);                     % Critical pressure
T_crit = prop_calculation('T_CRITICAL',working_fluid);                     % Critical temperature
% p_crit = [];
% T_crit = [];

% Define the triple properties (only pure substances)
p_trip = prop_calculation('P_TRIPLE',working_fluid);                       % Triple pressure
T_trip = prop_calculation('T_TRIPLE',working_fluid);                       % Triple temperature
% p_trip = [];
% T_trip = [];

% Define the saturation pressure and enthalpy at ambient temperature
% If the ambient temperature is higher than the critical temperature use
% the critical pressure and enthalpy instead
if T_0 < T_crit
    p_sat0 = prop_calculation('P','T',T_0,'Q',0,working_fluid);
    h_sat0 = prop_calculation('H','T',T_0,'Q',0,working_fluid);
else
    p_sat0 = p_crit;
    h_sat0 = prop_calculation('H','T',T_crit,'P',p_crit,working_fluid);
end

% Define the pressure and temperature limits
p_min = 1.01*p_trip;                                                       % Minimum pressure
p_max = 3.00*p_crit;                                                       % Maximum pressure
T_min = T_0;                                                               % Maximum temperature
T_max = 1.50*T_crit;                                                       % Maximum temperature

% Define the enthalpy limits
h_min = h_sat0;
h_max = prop_calculation('H','T',T_1,'P',1e-3,working_fluid);


%% Define the independent variables and bounds
% x1 = Heat source exit temperature
x1 = 0.25;
x1_min = 0.00;
x1_max = 1.00;

% x2 = Heat sink exit temperature
x2 = 0.25;
x2_min = 0.00; 
x2_max = 1.00;

% x3 = Pump inlet pressure (subcritical)
x3 = (1.5*p_sat0-p_min)/(p_max-p_min);
x3_min = 0.00;
x3_max = 1.00;

% x4 = Pump inlet enthalpy
x4 = 0.05;
x4_min = 0.00;
x4_max = 1.00;

% x5 = Expander inlet pressure
x5 = (0.50*p_crit-p_min)/(p_max-p_min);
x5_min = 0.00;
x5_max = 1.00;

% x6 = Expander inlet enthalpy
x6 = 0.80;
x6_min = 0.00;
x6_max = 1.00;

% x7 = Recuperator effectiveness (set to zero for no recuperator)
x7 = 0.00;
x7_min = 0.00;
x7_max = 0.00;

% Organize the information in vectors
lb_cycle = [x1_min, x2_min, x3_min, x4_min, x5_min, x6_min, x7_min]'; 
ub_cycle = [x1_max, x2_max, x3_max, x4_max, x5_max, x6_max, x7_max]';       
x0_cycle = [x1, x2, x3, x4, x5, x6, x7]';

% Some guidelines for the degrees of freedom
%{

Play around with x0_cycle to find a good initial guess.

Mathematical definition of the degrees of freedom
  x1 = (T_3 - T_3min) / (T_3max - T_3min)
  x2 = (T_12 - T_12min) / (T_12max - T_12min)
  x3 = (p_4 - p_min) / (p_max - p_min)
  x4 = (h_4 - h_min) / (h_max - h_min)
  x5 = (p_7 - p_min) / (p_max - p_min)
  x6 = (h_7 - h_min) / (h_max - h_min)
  x7 = (h6 - h5) / (h(p6,T8) - h5)

If you want a cycle configuration with no recuperator you have to:
 1) Set the recuperator effectiveness equal to zero:
    x7 = x7_min = x7_max = 0.00
 2) Set the pressure drops in the recuperator equal to zero:
    dp_h_rec = dp_c_rec = 0.00
 3) Do not apply a constraint for the recuperator pinch point dT:
    constraints.dT_rec.apply = 'no'

How to give start values when you use the code for the first time?
  x1 usually takes low values (maximum exploitation of the heat source)
  x2 usually takes low values (small temperature jump in the heat sink)
  x3 usually takes low values (low pressure at the exit of the condenser)
  x4 usually takes low values (low enthalpy at the exit of the condenser)
  x5<x_crit for subcritical cycles and x5>x_crit for transcritical cycles
  x6 usually takes high values (high enthalpy at the inlet of the evap.)
  x7=0 corresponds to no recuperation and x=1 to ideal recuperation

Note: Low means close to zero and high means close to one

%}


%% Define the constraints
% Instructions:
% 1. Specify the minimum and maximum values
% 2. Specify whether to apply the constraint or not
% 3. Use [brackets] to ignore the constraint value

% Minimum temperature difference in the evaporator (pinch point)
constraints.dT_evap.min = 10;
constraints.dT_evap.apply = 'yes';

% Minimum temperature difference in the condenser (pinch point)
constraints.dT_cond.min = 10;
constraints.dT_cond.apply = 'yes';

% Minimum temperature difference in the recuperator (pinch point)
constraints.dT_rec.min = 10;
constraints.dT_rec.apply = 'no';

% Minimum vapor quality along the expander
constraints.quality.min = 1.00;
constraints.quality.apply = 'yes';
 
% Degree of subcooling at the inlet of the pump
constraints.dT_subcooling.min = 1;
constraints.dT_subcooling.max = [];
constraints.dT_subcooling.apply = 'yes';

% Degree of superheating at inlet of the expander (subcritical cycles)
constraints.dT_superheating.min = 10;
constraints.dT_superheating.max = [];
constraints.dT_superheating.apply = 'yes';

% Maximum cycle pressure
constraints.pressure.min = [];
constraints.pressure.max = 20e5;
constraints.pressure.apply = 'yes';


%% Define optimization algorithm settings
% Set optimization algorithm
% 'spq' is my favourite because it converges reliably and respects bounds
% 'active-set' is a bit more aggressive and does not always respect the
% bounds. It is faster, but less reliable, than 'sqp'.
% 'interior-point' also works, but it requires more iterations
algorithm = 'sqp';

% Compute the problem gradients in parallel or not
use_parallel = false; % true or false

% Define termination criteria (no need to change in most cases)
max_iterations       = 1000;
max_function_evals   = 10000;
step_tolerance       = 1e-08;
function_tolerance   = 1e-05;
constraint_tolerance = 1e-05;
optimality_tolerance = 1e-05;   
                                  
% Description of the possible outputs of the optimization function
%{

'interior-point', 'active-set', and 'sqp' algorithms
Exitflag =  1 -> Success. First-order optimality measure was less than options and constraints are not violated.
Exitflag =  2 -> Success. Step size was less than options and constraints are not violated.
Exitflag =  0 -> Unsuccess. Number of iterations or function evaluations exceeded option
Exitflag = -1 -> Unsuccess. Solver stopped by an output function
Exitflag = -2 -> Unsuccess. No feasible point was found

Only active set algorithm
Exitflag =  4 -> Success. Magnitude of the search direction was less than 2*options.StepTolerance and maximum constraint violation was less than options.ConstraintTolerance.
Exitflag =  5 -> Success. Magnitude of directional derivative in search direction was less than 2*options.OptimalityTolerance and maximum constraint violation was less than options.ConstraintTolerance.

%}


%% Choose what figures should be plotted
choose_plots = struct('diagram_TQ',            'yes', ...                  % Temperature - heat duty diagram
                      'diagram_Ts',            'yes', ...                  % Temperature - entropy diagram of the cycle
                      'diagram_Th',            'yes', ...                   % Temperature - enthalpy diagram of the cycle
                      'save',                  0);                         % Choose whether to save the figures or not
                  
% Description of the saving options
%{

save=0 does not save the figures
save=1 saves the figure in vector format using the default MATLAB export function (faster)
save=2 saves the figure in vector using the export_fig library (requires ghostscript)     
                                                              
In order to use save=2 it is necessary to install 'ghoscript'
  https://github.com/altmany/export_fig
  http://www.ghostscript.com
  http://xpdfreader.com

%}


%% Store the parameters defined up to this point into data structures
% The script 'create_problem_structures' located in source_code folder
%  1) 'fixed_parameters' stores the cycle specifications 
%  2) 'optimization_problem' stores the optimization problem settings
create_problem_structures

% These data structures are created  in a different script to keep the
% number of lines of this script to a minimum

% Be tidy and clear the worskpace now that we stored everything
clearvars -except optimization_problem fixed_parameters project_name results_path
                   

%% Load a previous solution as initial guess
% % Uncomment these lines to use a previous solution as initial guess
% load(fullfile(fixed_parameters.results_path,'cycle_data.mat'))
% optimization_problem.x0 = cycle_data.optimization.x;


%% Plot the initial solution
filename_suffix = 'initial';
cycle_data = create_plots(optimization_problem.x0,fixed_parameters,filename_suffix);
print_solution(cycle_data)


%% Solve the optimization problem
[cycle_data,x_opt,f_opt,exitflag,output,lambda] = solve_optimization_problem(fixed_parameters,optimization_problem);
save_solution(cycle_data, project_name, results_path)

%% Plot the optimal solution
filename_suffix = 'optimal';
create_plots(x_opt,fixed_parameters,filename_suffix);
print_solution(cycle_data)



