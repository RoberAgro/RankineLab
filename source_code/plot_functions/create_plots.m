function cycle_data = create_plots(x,fixed_parameters,filename_suffix)

% Load setting for beautiful plots
set_plot_options()

% Evaluate the optimization problem
fixed_parameters.calc_detail = 'long';
if nargin == 2
    filename = fixed_parameters.project_name;
elseif nargin == 3
    filename = [fixed_parameters.project_name, '_', filename_suffix];
else
    error('The number of arguments must be 2 or 3')
end
filepath = fixed_parameters.results_path;
cycle_data = evaluate_optimization_problem(x,fixed_parameters);

% Choose whether to save the plots or not
% 'choose_plots' is an structure that contains what diagrams to draw
save = fixed_parameters.choose_plots.save;

% Plot temperature vs heat flow rate diagram
if strcmp(fixed_parameters.choose_plots.diagram_TQ,'yes')
    plot_TQ_diagram(cycle_data,filename,filepath,save)
end

% Plot temperature vs entropy diagram
if strcmp(fixed_parameters.choose_plots.diagram_Ts,'yes')
    plot_Ts_diagram(cycle_data,filename,filepath,save)
end

% Plot temperature vs enthalpy diagram
if strcmp(fixed_parameters.choose_plots.diagram_Th,'yes')
    plot_Th_diagram(cycle_data,filename,filepath,save)
end

end
