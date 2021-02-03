function [cycle_data, f, c, c_eq] = evaluate_optimization_problem(x,fixed_parameters)

% Evaluate the thermodynamic cycle model
cycle_data = evaluate_rankine_cycle(x,fixed_parameters);

% Extract the objective function value
f = -[cycle_data.efficiency.eta_1];
% f = -[cycle_data.efficiency.eta_2];

% Extract the vector of inequality constraints
c = cycle_data.optimization.c;

% Extract the vector of  equality constraints
c_eq = cycle_data.optimization.c_eq;
 
end

