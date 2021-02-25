%% Save the cycle parameters as a structure
% Simple script to organize information into an structure

% Thermodynamic boundaries
boundary_conditions = struct('heating_fluid',   {heating_fluid},   ...
                             'working_fluid',   {working_fluid},   ...
                             'cooling_fluid',   {cooling_fluid},   ...
                             'm_h',             m_h,               ...
                             'T_1',             T_1,               ...
                             'p_1',             p_1,               ...
                             'p_3',             p_3,               ...
                             'T_3min',          T_3min,            ...
                             'T_3max',          T_3max,            ...
                             'T_10',            T_10,              ...
                             'p_10',            p_10,              ...
                             'p_12',            p_12,              ...
                             'T_12min',         T_12min,            ...
                             'T_12max',         T_12max);

% Fluid properties
properties = struct('p_0',      p_0,      ...
                    'T_0',      T_0,      ...
                    'p_crit',   p_crit,   ...
                    'T_crit',   T_crit,   ...
                    'p_trip',   p_trip,   ...
                    'T_trip',   T_trip,   ...
                    'p_max',    p_max,    ...
                    'T_max',    T_max,    ...
                    'h_max',    h_max,    ...
                    'p_min',    p_min,    ...
                    'T_min',    T_min,    ...
                    'h_min',    h_min);

% Component specifications                     
components = struct('dp_h_evap',    dp_h_evap,     ...
                    'dp_c_evap',    dp_c_evap,     ...
                    'dp_h_cond',    dp_h_cond,     ...
                    'dp_c_cond',    dp_c_cond,     ...
                    'dp_h_rec',     dp_h_rec,      ...
                    'dp_c_rec',     dp_c_rec,      ...
                    'N_evap',       N_evap,        ...
                    'N_cond',       N_cond,        ...
                    'N_rec',        N_rec,         ...
                    'eta_pump_h',   eta_pump_h,     ...
                    'eta_pump_f',   eta_pump_f,     ...
                    'eta_pump_c',   eta_pump_c,     ...
                    'eta_expander', eta_expander,  ...
                    'pump_efficiency_definition', pump_efficiency_definition, ...
                    'expander_efficiency_definition', expander_efficiency_definition);
  
% Lower bounds for the degrees of freedom                 
lower_bounds = struct('x1',   x1_min,   ...
                      'x2',   x2_min,   ...
                      'x3',   x3_min,   ...
                      'x4',   x4_min,   ...
                      'x5',   x5_min,   ...
                      'x6',   x6_min,   ...
                      'x7',   x7_min);

% Upper bounds for the degrees of freedom     
upper_bounds = struct('x1',   x1_max,   ...
                      'x2',   x2_max,   ...
                      'x3',   x3_max,   ...
                      'x4',   x4_max,   ...
                      'x5',   x5_max,   ...
                      'x6',   x6_max,   ...
                      'x7',   x7_max);

% Organize the previous structures in a compact way     
fixed_parameters = struct('boundary_conditions', boundary_conditions,   ...
                          'properties',          properties,            ...
                          'components',          components,            ...
                          'lower_bounds',        lower_bounds,          ...
                          'upper_bounds',        upper_bounds,          ...
                          'constraints',         constraints,           ...
                          'results_path',        results_path,          ...
                          'project_name',        project_name,          ...
                          'choose_plots',        choose_plots);
                    
                      
%% Create optimization_problem structure
optimization_problem = struct('lb',        [], ...                         % Vector of lower bounds
                              'ub',        [], ...                         % Vector of upper bounds
                              'Aeq',       [], ...                         % Matrix for linear equality constraints
                              'beq',       [], ...                         % Vector for linear equality constraints
                              'Aineq',     [], ...                         % Matrix for linear inequality constraints
                              'bineq',     [], ...                         % Vector for linear inequality constraints
                              'nonlcon',   [], ...                         % Nonlinear constraints function
                              'objective', [], ...                         % Objective function
                              'x0',        [], ...                         % Initial guess for the degrees of freedom
                              'options',   [], ...                         % Optimization options structure
                              'solver', 'fmincon');

% Define the initial guess vector
optimization_problem.x0 = x0_cycle;

% Define the vector of lower bounds
optimization_problem.lb = lb_cycle;
                       
% Define the vector of upper bounds                  
optimization_problem.ub = ub_cycle;

% Define the options using the optimoptions function
% Change the type of algorithm here
optimization_problem.options = optimoptions(@fmincon,                   ...
                       'Display', 'iter-detailed',                      ...
                       'Algorithm', algorithm,                          ...
                       'StepTolerance', step_tolerance,                 ...
                       'ConstraintTolerance', constraint_tolerance,     ...
                       'OptimalityTolerance', optimality_tolerance,     ...
                       'FunctionTolerance', function_tolerance,         ...
                       'MaxFunctionEvaluations', max_function_evals,    ...
                       'MaxIterations', max_iterations,                 ...
                       'FiniteDifferenceStepSize', sqrt(eps),           ...
                       'UseParallel', use_parallel,                     ...
                       'PlotFcn', {@optimplotfval,                      ...
                                   @optimplotconstrviolation,           ...
                                   @optimplotfirstorderopt,             ...
                                   @optimplotstepsize});
              

optimization_problem.options.OutputFcn = @(x,optimValues,state)save_current_solution(x,optimValues,state,fixed_parameters);
    


