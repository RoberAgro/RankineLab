function cycle_data = evaluate_rankine_cycle(x,fixed_parameters)

% This function contains the mathematical model of the thermodynamic cycle

%% Load parameters
% Load working fluids
heating_fluid = fixed_parameters.boundary_conditions.heating_fluid;
working_fluid = fixed_parameters.boundary_conditions.working_fluid;
cooling_fluid = fixed_parameters.boundary_conditions.cooling_fluid;

% Load heat source specifications
m_h = fixed_parameters.boundary_conditions.m_h;
T_1 = fixed_parameters.boundary_conditions.T_1;
p_1 = fixed_parameters.boundary_conditions.p_1;
p_3 = fixed_parameters.boundary_conditions.p_3;
T_3max = fixed_parameters.boundary_conditions.T_3max;
T_3min = fixed_parameters.boundary_conditions.T_3min;

% Load the heat sink specifications
T_10 = fixed_parameters.boundary_conditions.T_10;
p_10 = fixed_parameters.boundary_conditions.p_10;
p_12 = fixed_parameters.boundary_conditions.p_12;
T_12max = fixed_parameters.boundary_conditions.T_12max;
T_12min = fixed_parameters.boundary_conditions.T_12min;

% Load pressure drops in the heat exchangers
dp_h_evap = fixed_parameters.components.dp_h_evap;
dp_c_evap = fixed_parameters.components.dp_c_evap;
dp_h_cond = fixed_parameters.components.dp_h_cond;
dp_c_cond = fixed_parameters.components.dp_c_cond;
dp_h_rec = fixed_parameters.components.dp_h_rec;
dp_c_rec = fixed_parameters.components.dp_c_rec;

% Load the number of discretizations of the heat exchangers
N_evap = fixed_parameters.components.N_evap;
N_cond = fixed_parameters.components.N_cond;
N_rec = fixed_parameters.components.N_rec;

% Load turbomachinery efficiency
eta_pump_h = fixed_parameters.components.eta_pump_h;
eta_pump_f = fixed_parameters.components.eta_pump_f;
eta_pump_c = fixed_parameters.components.eta_pump_c;
eta_expander = fixed_parameters.components.eta_expander;
eta_expander_definition = fixed_parameters.components.eta_expander_definition;

% Load fluid properties
p_0 = fixed_parameters.properties.p_0;
T_0 = fixed_parameters.properties.T_0;
p_crit = fixed_parameters.properties.p_crit;
T_crit = fixed_parameters.properties.T_crit;
p_trip = fixed_parameters.properties.p_trip;
T_trip = fixed_parameters.properties.T_trip;
p_max = fixed_parameters.properties.p_max;
T_max = fixed_parameters.properties.T_max;
h_max = fixed_parameters.properties.h_max;
p_min = fixed_parameters.properties.p_min;
T_min = fixed_parameters.properties.T_min;
h_min = fixed_parameters.properties.h_min;

% Load the level of detail in the coumputations
% 'short' for optimization
% 'long' for plotting
calc_detail= fixed_parameters.calc_detail;


%% Degrees of freedom
T_3  = T_3min+x(1)*(T_3max-T_3min);
T_12 = T_12min+x(2)*(T_12max-T_12min);
p_4  = p_min+x(3)*(p_max-p_min);
h_4  = h_min+x(4)*(h_max-h_min);
p_7  = p_min+x(5)*(p_max-p_min);
h_7  = h_min+x(6)*(h_max-h_min);
effectiveness = x(7);

% % Nasty trick to avoid computations too close to the critical point
% if p_4/p_crit > 0.99 && p_4/p_crit < 1.01
%     p_4 = 0.99*p_crit;
%     disp('Warning: p_4 too close to the critical point, p_4 set to 0.99*p_crit')
% end
% 
% if p_7/p_crit > 0.99 && p_7/p_crit < 1.01
%     p_7 = 0.99*p_crit;
%     disp('Warning: p_7 too close to the critical point, p_7 set to 0.99*p_crit')
% end

%% Sequential computation of the states
% State 1 (heat source inlet)
h_1 = prop_calculation('H','T',T_1,'P',p_1,heating_fluid);
d_1 = prop_calculation('D','T',T_1,'P',p_1,heating_fluid);
s_1 = prop_calculation('S','T',T_1,'P',p_1,heating_fluid);

% State 3 (heat source outlet)
h_3 = prop_calculation('H','T',T_3,'P',p_3,heating_fluid);
d_3 = prop_calculation('D','T',T_3,'P',p_3,heating_fluid);
s_3 = prop_calculation('S','T',T_3,'P',p_3,heating_fluid);

% State 2 (heat source pump inlet)
p_2 = p_1*(1-dp_h_evap);
[h_2,p_pump_h,h_pump_h,T_pump_h,s_pump_h,d_pump_h] = evaluate_compressor(heating_fluid,h_3,p_3,p_2,eta_pump_h,calc_detail);
T_2 = prop_calculation('T','H',h_2,'P',p_2,heating_fluid);
s_2 = prop_calculation('S','H',h_2,'P',p_2,heating_fluid);
d_2 = prop_calculation('D','H',h_2,'P',p_2,heating_fluid);

% State 3 (pump inlet)
T_4 = prop_calculation('T','P',p_4,'H',h_4,working_fluid);
s_4 = prop_calculation('S','P',p_4,'H',h_4,working_fluid);
d_4 = prop_calculation('D','P',p_4,'H',h_4,working_fluid);

% State 5 (pump outlet)
p_6 = p_7/(1-dp_c_evap);
p_5 = p_6/(1-dp_c_rec);
[h_5,p_pump_f,h_pump_f,T_pump_f,s_pump_f,d_pump_f] = evaluate_compressor(working_fluid,h_4,p_4,p_5,eta_pump_f,calc_detail);
T_5 = prop_calculation('T','P',p_5,'H',h_5,working_fluid);
s_5 = prop_calculation('S','P',p_5,'H',h_5,working_fluid);
d_5 = prop_calculation('D','P',p_5,'H',h_5,working_fluid);

% State 7 (expander inlet)
T_7 = prop_calculation('T','P',p_7,'H',h_7,working_fluid);
s_7 = prop_calculation('S','P',p_7,'H',h_7,working_fluid);
d_7 = prop_calculation('D','P',p_7,'H',h_7,working_fluid);

% State 8 (expander outlet)
p_9 = p_4/(1-dp_h_cond);
p_8 = p_9/(1-dp_h_rec);
[h_8,p_exp,h_exp,T_exp,s_exp,d_exp] = evaluate_expander(working_fluid,h_7,p_7,p_8,eta_expander,eta_expander_definition,calc_detail);
T_8 = prop_calculation('T','P',p_8,'H',h_8,working_fluid);
s_8 = prop_calculation('S','P',p_8,'H',h_8,working_fluid);
d_8 = prop_calculation('D','P',p_8,'H',h_8,working_fluid);

% State 5 (evaporator inlet)
try
    h_6ideal = prop_calculation('H','T',T_8,'P',p_6,working_fluid);
catch
    h_6ideal = prop_calculation('H','T',T_8,'Q',0,working_fluid);
    disp('Warning: h_8ideal computation failed, using h_8sat instead');
end
h_6 = h_5 + effectiveness * (h_6ideal - h_5); % Effectiveness definition
T_6 = prop_calculation('T','P',p_6,'H',h_6,working_fluid);
s_6 = prop_calculation('S','P',p_6,'H',h_6,working_fluid);
d_6 = prop_calculation('D','P',p_6,'H',h_6,working_fluid);

% State 9 (condenser inlet)
h_9 = h_8 - (h_6 - h_5); % Energy balance recuperator
T_9 = prop_calculation('T','P',p_9,'H',h_9,working_fluid);
s_9 = prop_calculation('S','P',p_9,'H',h_9,working_fluid);
d_9 = prop_calculation('D','P',p_9,'H',h_9,working_fluid);

% State 10 (heat sink inlet)
h_10 = prop_calculation('H','T',T_10,'P',p_10,cooling_fluid);
s_10 = prop_calculation('S','T',T_10,'P',p_10,cooling_fluid);
d_10 = prop_calculation('D','T',T_10,'P',p_10,cooling_fluid);

% State 11 (heat sink pump outlet)
p_11 = p_12/(1-dp_c_cond);
[h_11,p_pump_c,h_pump_c,T_pump_c,s_pump_c,d_pump_c] = evaluate_compressor(cooling_fluid,h_10,p_10,p_11,eta_pump_c,calc_detail);
T_11 = prop_calculation('T','H',h_11,'P',p_11,cooling_fluid);
s_11 = prop_calculation('S','H',h_11,'P',p_11,cooling_fluid);
d_11 = prop_calculation('D','H',h_11,'P',p_11,cooling_fluid);

% State 11 (heat sink outlet)
h_12 = prop_calculation('H','T',T_12,'P',p_12,cooling_fluid);
s_12 = prop_calculation('S','T',T_12,'P',p_12,cooling_fluid);
d_12 = prop_calculation('D','T',T_12,'P',p_12,cooling_fluid);


%% Check the quality of the fluid along the expansion
q_exp = 0*p_exp;
for i = 1:length(p_exp)
    q_exp(i) = quality('H',h_exp(i),working_fluid,p_exp(i),p_crit);
end


%% Compute the mass flow rates
m_f = m_h*(h_1-h_2)/(h_7-h_6);        % Working fluid mass flow rate
m_c = m_f*(h_9-h_4)/(h_12-h_11);      % Cooling fluid mass flow rate


%% Compute the pinch point in the heat exchangers
[dT_evap,T_evap_h,T_evap_c, ...
         h_evap_h,h_evap_c, ...
         p_evap_h,p_evap_c, ...
         d_evap_h,d_evap_c, ...
         s_evap_h,s_evap_c] = evaluate_exchanger(heating_fluid,working_fluid,m_h,m_f,p_2,p_1,p_6,p_7,h_2,h_1,h_6,h_7,N_evap,calc_detail);

[dT_cond,T_cond_h,T_cond_c, ...
         h_cond_h,h_cond_c, ...
         p_cond_h,p_cond_c, ...
         d_cond_h,d_cond_c, ...
         s_cond_h,s_cond_c] = evaluate_exchanger(working_fluid,cooling_fluid,m_f,m_c,p_4,p_9,p_11,p_12,h_4,h_9,h_11,h_12,N_cond,calc_detail);

[dT_rec,T_rec_h,T_rec_c, ...
         h_rec_h,h_rec_c, ...
         p_rec_h,p_rec_c, ...
         d_rec_h,d_rec_c, ...
         s_rec_h,s_rec_c] = evaluate_exchanger(working_fluid,working_fluid,m_f,m_f,p_9,p_8,p_5,p_6,h_9,h_8,h_5,h_6,N_rec,calc_detail);

        
%%  Dead states
h_h0 = prop_calculation('H','T',T_0,'P',p_0,heating_fluid);
s_h0 = prop_calculation('S','T',T_0,'P',p_0,heating_fluid);
h_f0 = prop_calculation('H','T',T_0,'P',p_0,working_fluid);
s_f0 = prop_calculation('S','T',T_0,'P',p_0,working_fluid);
h_c0 = prop_calculation('H','T',T_0,'P',p_0,cooling_fluid);
s_c0 = prop_calculation('S','T',T_0,'P',p_0,cooling_fluid);


%% Energy analysis
% Heat and power flows
Q_max  = m_h*(h_1-h_h0);
Q_evap = m_h*(h_1-h_2);
Q_cond = m_f*(h_9-h_4);
Q_rec = m_f*(h_8-h_9);
W_exp = m_f*(h_7-h_8);
W_pump_h = m_h*(h_3-h_2);
W_pump_f = m_f*(h_5-h_4);
W_pump_c = m_c*(h_11-h_10);
W_net = W_exp-W_pump_h-W_pump_f-W_pump_c;

% Flow energy in an out of the power cycle
E_h_in = m_h*(h_1-h_h0);
E_h_out = m_h*(h_3-h_h0);
E_c_in = m_c*(h_10-h_c0);
E_c_out = m_c*(h_12-h_c0);

% Check if the energy balance matches
energy_check1 = [E_h_in, -E_h_out, E_c_in, -E_c_out, W_pump_h, W_pump_f, W_pump_c, -W_exp] / E_h_in;
energy_check2 = [Q_evap, W_pump_f, -W_exp, -Q_cond]/E_h_in;
energy_error1 = sum(energy_check1);
energy_error2 = sum(energy_check2);

% Efficiency definitions
eta_1recovery = Q_evap/Q_max;
eta_1cycle    = W_net/Q_evap;
eta_1plant    = W_net/Q_max;


%% Exergy analysis
% Heat source exergy flows
e_1 = (h_1-h_h0)-T_0*(s_1-s_h0);
e_2 = (h_2-h_h0)-T_0*(s_2-s_h0);
e_3 = (h_3-h_h0)-T_0*(s_3-s_h0);
e_4 = (h_4-h_f0)-T_0*(s_4-s_f0);
e_5 = (h_5-h_f0)-T_0*(s_5-s_f0);
e_6 = (h_6-h_f0)-T_0*(s_6-s_f0);
e_7 = (h_7-h_f0)-T_0*(s_7-s_f0);
e_8 = (h_8-h_f0)-T_0*(s_8-s_f0);
e_9 = (h_9-h_f0)-T_0*(s_9-s_f0);
e_10 = (h_10-h_c0)-T_0*(s_10-s_c0);
e_11 = (h_11-h_c0)-T_0*(s_11-s_c0);
e_12 = (h_12-h_c0)-T_0*(s_12-s_c0);

% Exergy destruction in each component
I_pump_h = m_h*T_0*(s_3-s_2);
I_pump_f = m_f*T_0*(s_5-s_4);
I_pump_c = m_c*T_0*(s_11-s_10);
I_exp = m_f*T_0*(s_8-s_7);
I_evap = m_h*(e_1-e_2)-m_f*(e_7-e_6);
I_cond = m_f*(e_9-e_4)-m_c*(e_12-e_11);
I_rec = m_f*(e_8-e_9)-m_f*(e_6-e_5);

% Flow exergy in an out of the power cycle
Ex_h_in = m_h*e_1;
Ex_h_out = m_h*e_3;
Ex_c_in = m_c*e_10;
Ex_c_out = m_c*e_12;
Ex_evap = m_h*(e_1-e_2);

% Efficiency definitions
eta_2recovery = Ex_evap/Ex_h_in;
eta_2cycle    = W_net/Ex_evap;
eta_2plant    = W_net/Ex_h_in;
eta_1max      = Ex_h_in/E_h_in;

% Check exergy computations
exergy_check = [Ex_h_in, -Ex_h_out, Ex_c_in, -Ex_c_out, W_pump_h, W_pump_f, W_pump_c, -W_exp, -I_exp, -I_pump_h, -I_pump_f, -I_pump_c, -I_evap, -I_cond, -I_rec]/Ex_h_in;
exergy_error = sum(exergy_check);


%% Evaluate the optimization problem constraints     
% Initialize constraint vectors
c_eq = [];      % Vector of equality constraints
c_ineq = [];    % Vector of inequality constraints

% Limitation for the minimum temperature difference along the evaporator
dT_evap_min = fixed_parameters.constraints.dT_evap.min;  
if strcmp(fixed_parameters.constraints.dT_evap.apply,'yes') == 1
    c_ineq = [c_ineq; -(dT_evap-dT_evap_min)];
end

% Limitation for the minimum temperature difference along the condenser
dT_cond_min = fixed_parameters.constraints.dT_cond.min;
if strcmp(fixed_parameters.constraints.dT_cond.apply,'yes') == 1
    c_ineq = [c_ineq; -(dT_cond-dT_cond_min)];
end

% Limitation for the minimum temperature difference along the recuperator
dT_rec_min = fixed_parameters.constraints.dT_rec.min;  
if strcmp(fixed_parameters.constraints.dT_rec.apply,'yes') == 1
    c_ineq = [c_ineq; -(dT_rec-dT_rec_min)];
end

% Limitation for the vapor quality along the expander
expander_quality_min = fixed_parameters.constraints.quality.min;
if strcmp(fixed_parameters.constraints.quality.apply,'yes') == 1
    c_ineq = [c_ineq; -(q_exp-expander_quality_min)];
end

% Limitation for the min/max subcooling at the inlet of the pump
% It is necessary to compute it in terms of enthalpy changes because the
% temperature would be insensitive to changes in the degrees of freedom if
% the fluid is in the two-phase region (the optimizer gets stuck)
dT_subcooling_min = fixed_parameters.constraints.dT_subcooling.min;
dT_subcooling_max = fixed_parameters.constraints.dT_subcooling.max;
if p_4 < p_crit
    T_4sat = prop_calculation('T','P',p_4,'Q',0,working_fluid);
else
    T_4sat = T_crit;
end
dT_subcooling = T_4sat - T_4;
if strcmp(fixed_parameters.constraints.dT_subcooling.apply,'yes') == 1
    if isempty(dT_subcooling_min) ~= 1
        h_4_min = prop_calculation('H','T',T_4sat-dT_subcooling_min,'P',p_4,working_fluid);
        c_ineq = [c_ineq; +(h_4-h_4_min)/(h_max-h_min)];    % Careful with sign, confusing
    end
    if isempty(dT_subcooling_max) ~= 1
        h_4_max = prop_calculation('H','T',T_4sat-dT_subcooling_max,'P',p_4,working_fluid);
        c_ineq = [c_ineq; -(h_4-h_4_max)/(h_max-h_min)];    % Careful with sign, confusing
    end
end

% Limitation for the min/max superheating at the inlet of the expander
% It is necessary to compute it in terms of enthalpy changes because the
% temperature would be insensitive to changes in the degrees of freedom if
% the fluid is in the two-phase region (the optimizer gets stuck)
% Limitation for the min/max superheating at the inlet of the expander
dT_superheating_min = fixed_parameters.constraints.dT_superheating.min;
dT_superheating_max = fixed_parameters.constraints.dT_superheating.max;
if p_7 < p_crit
    T_7sat = prop_calculation('T','P',p_7,'Q',1,working_fluid);
else
    T_7sat = T_crit;
end
dT_superheating = T_7 - T_7sat;
if strcmp(fixed_parameters.constraints.dT_superheating.apply,'yes') == 1
    if isempty(dT_superheating_min) ~= 1
        h_7_min = prop_calculation('H','T',T_7sat+dT_superheating_min,'P',p_7,working_fluid);
        c_ineq = [c_ineq; -(h_7-h_7_min)/(h_max-h_min)];
    end
    if isempty(dT_superheating_max) ~= 1
        h_7_max = prop_calculation('H','T',T_7sat+dT_superheating_max,'P',p_7,working_fluid);
        c_ineq = [c_ineq; +(h_7-h_7_max)/(h_max-h_min)];
    end
end

% Limitation for the min/max pressure of the working fluid
p_cycle = [p_4, p_5, p_6, p_7, p_8, p_9]';
p_min_constraint = fixed_parameters.constraints.pressure.min;
p_max_constraint = fixed_parameters.constraints.pressure.max;
if strcmp(fixed_parameters.constraints.pressure.apply,'yes') == 1
    if isempty(p_min_constraint) ~= 1
        c_ineq = [c_ineq; -(p_cycle-p_min_constraint)];
    end
    if isempty(p_max_constraint) ~= 1
        c_ineq = [c_ineq; +(p_cycle-p_max_constraint)];
    end
end


%% Export the cycle variables
% Create the data structure
cycle_data = create_data_structure;    % Function without inputs

% Temperatures
cycle_data.cycle_states(1).fluid = heating_fluid;
cycle_data.cycle_states(2).fluid = heating_fluid;
cycle_data.cycle_states(3).fluid = heating_fluid;
cycle_data.cycle_states(4).fluid = working_fluid;
cycle_data.cycle_states(5).fluid = working_fluid;
cycle_data.cycle_states(6).fluid = working_fluid;
cycle_data.cycle_states(7).fluid = working_fluid;
cycle_data.cycle_states(8).fluid = working_fluid;
cycle_data.cycle_states(9).fluid = working_fluid;
cycle_data.cycle_states(10).fluid = cooling_fluid;
cycle_data.cycle_states(11).fluid = cooling_fluid;
cycle_data.cycle_states(12).fluid = cooling_fluid;

% Temperatures
cycle_data.cycle_states(1).T = T_1;
cycle_data.cycle_states(2).T = T_2;
cycle_data.cycle_states(3).T = T_3;
cycle_data.cycle_states(4).T = T_4;
cycle_data.cycle_states(5).T = T_5;
cycle_data.cycle_states(6).T = T_6;
cycle_data.cycle_states(7).T = T_7;
cycle_data.cycle_states(8).T = T_8;
cycle_data.cycle_states(9).T = T_9;
cycle_data.cycle_states(10).T = T_10;
cycle_data.cycle_states(11).T = T_11;
cycle_data.cycle_states(12).T = T_12;

% Presures
cycle_data.cycle_states(1).p = p_1;
cycle_data.cycle_states(2).p = p_2;
cycle_data.cycle_states(3).p = p_3;
cycle_data.cycle_states(4).p = p_4;
cycle_data.cycle_states(5).p = p_5;
cycle_data.cycle_states(6).p = p_6;
cycle_data.cycle_states(7).p = p_7;
cycle_data.cycle_states(8).p = p_8;
cycle_data.cycle_states(9).p = p_9;
cycle_data.cycle_states(10).p = p_10;
cycle_data.cycle_states(11).p = p_11;
cycle_data.cycle_states(12).p = p_12;

% Densities
cycle_data.cycle_states(1).d = d_1;
cycle_data.cycle_states(2).d = d_2;
cycle_data.cycle_states(3).d = d_3;
cycle_data.cycle_states(4).d = d_4;
cycle_data.cycle_states(5).d = d_5;
cycle_data.cycle_states(6).d = d_6;
cycle_data.cycle_states(7).d = d_7;
cycle_data.cycle_states(8).d = d_8;
cycle_data.cycle_states(9).d = d_9;
cycle_data.cycle_states(10).d = d_10;
cycle_data.cycle_states(11).d = d_11;
cycle_data.cycle_states(12).d = d_12;

% Enthalpies
cycle_data.cycle_states(1).h = h_1;
cycle_data.cycle_states(2).h = h_2;
cycle_data.cycle_states(3).h = h_3;
cycle_data.cycle_states(4).h = h_4;
cycle_data.cycle_states(5).h = h_5;
cycle_data.cycle_states(6).h = h_6;
cycle_data.cycle_states(7).h = h_7;
cycle_data.cycle_states(8).h = h_8;
cycle_data.cycle_states(9).h = h_9;
cycle_data.cycle_states(10).h = h_10;
cycle_data.cycle_states(11).h = h_11;
cycle_data.cycle_states(12).h = h_12;

% Entropies
cycle_data.cycle_states(1).s = s_1;
cycle_data.cycle_states(2).s = s_2;
cycle_data.cycle_states(3).s = s_3;
cycle_data.cycle_states(4).s = s_4;
cycle_data.cycle_states(5).s = s_5;
cycle_data.cycle_states(6).s = s_6;
cycle_data.cycle_states(7).s = s_7;
cycle_data.cycle_states(8).s = s_8;
cycle_data.cycle_states(9).s = s_9;
cycle_data.cycle_states(10).s = s_10;
cycle_data.cycle_states(11).s = s_11;
cycle_data.cycle_states(12).s = s_12;

% Exergies
cycle_data.cycle_states(1).e = e_1;
cycle_data.cycle_states(2).e = e_2;
cycle_data.cycle_states(3).e = e_3;
cycle_data.cycle_states(4).e = e_4;
cycle_data.cycle_states(5).e = e_5;
cycle_data.cycle_states(6).e = e_6;
cycle_data.cycle_states(7).e = e_7;
cycle_data.cycle_states(8).e = e_8;
cycle_data.cycle_states(9).e = e_9;
cycle_data.cycle_states(10).e = e_10;
cycle_data.cycle_states(11).e = e_11;
cycle_data.cycle_states(12).e = e_12;

% Mass flow rates
cycle_data.mass_flows.m_h = m_h;
cycle_data.mass_flows.m_f = m_f;
cycle_data.mass_flows.m_c = m_c;

% Fluids
cycle_data.fluids.heating_fluid = heating_fluid;
cycle_data.fluids.working_fluid = working_fluid;
cycle_data.fluids.cooling_fluid = cooling_fluid;

% Energy flows
cycle_data.energy_flows.Q_max  = Q_max;
cycle_data.energy_flows.Q_evap = Q_evap;
cycle_data.energy_flows.Q_cond = Q_cond;
cycle_data.energy_flows.Q_rec = Q_rec;
cycle_data.energy_flows.W_exp = W_exp;
cycle_data.energy_flows.W_pump_h = W_pump_h;
cycle_data.energy_flows.W_pump_f = W_pump_f;
cycle_data.energy_flows.W_pump_c = W_pump_c;

% Energy fractions
cycle_data.energy_frac.E_h_in    = 1.00;
cycle_data.energy_frac.E_h_out   = E_h_out/E_h_in;
cycle_data.energy_frac.E_c_in    = E_c_in/E_h_in;
cycle_data.energy_frac.E_c_out   = E_c_out/E_h_in;
cycle_data.energy_frac.W_exp     = W_exp/E_h_in;
cycle_data.energy_frac.W_pump_h  = W_pump_h/E_h_in;
cycle_data.energy_frac.W_pump_f  = W_pump_f/E_h_in;
cycle_data.energy_frac.W_pump_c  = W_pump_c/E_h_in;
cycle_data.energy_frac.Q_evap    = Q_evap/E_h_in;
cycle_data.energy_frac.Q_cond    = Q_cond/E_h_in;
cycle_data.energy_frac.Q_rec     = Q_rec/E_h_in;
cycle_data.energy_frac.W_net     = W_net/E_h_in;

% Exergy flows
cycle_data.exergy_flows.Ex_h_in   = Ex_h_in;
cycle_data.exergy_flows.Ex_h_out  = Ex_h_out;
cycle_data.exergy_flows.Ex_c_in   = Ex_c_in;
cycle_data.exergy_flows.Ex_c_out  = Ex_c_out;
cycle_data.exergy_flows.I_exp     = I_exp;
cycle_data.exergy_flows.I_pump_h  = I_pump_h;
cycle_data.exergy_flows.I_pump_f  = I_pump_f;
cycle_data.exergy_flows.I_pump_c  = I_pump_c;
cycle_data.exergy_flows.I_evap    = I_evap;
cycle_data.exergy_flows.I_cond    = I_cond;
cycle_data.exergy_flows.I_rec     = I_rec;

% Exergy fractions
cycle_data.exergy_frac.Ex_h_in   = 1.00;
cycle_data.exergy_frac.Ex_h_out  = Ex_h_out/Ex_h_in;
cycle_data.exergy_frac.Ex_c_in   = Ex_c_in/Ex_h_in;
cycle_data.exergy_frac.Ex_c_out  = Ex_c_out/Ex_h_in;
cycle_data.exergy_frac.I_exp     = I_exp/Ex_h_in;
cycle_data.exergy_frac.I_pump_h  = I_pump_h/Ex_h_in;
cycle_data.exergy_frac.I_pump_f  = I_pump_f/Ex_h_in;
cycle_data.exergy_frac.I_pump_c  = I_pump_c/Ex_h_in;
cycle_data.exergy_frac.I_evap    = I_evap/Ex_h_in;
cycle_data.exergy_frac.I_cond    = I_cond/Ex_h_in;
cycle_data.exergy_frac.I_rec     = I_rec/Ex_h_in;

% Cycle efficiencies
cycle_data.efficiency.eta_1 = eta_1plant;
cycle_data.efficiency.eta_2 = eta_2plant;
cycle_data.efficiency.eta_1recovery = eta_1recovery;
cycle_data.efficiency.eta_2recovery = eta_2recovery;
cycle_data.efficiency.eta_1cycle = eta_1cycle;
cycle_data.efficiency.eta_2cycle = eta_2cycle;
cycle_data.efficiency.eta_1max = eta_1max;

% Expander thermodynamic trajectory
cycle_data.expander.p = p_exp;
cycle_data.expander.h = h_exp;
cycle_data.expander.T = T_exp;
cycle_data.expander.s = s_exp;
cycle_data.expander.d = d_exp;
cycle_data.expander.q = q_exp;
cycle_data.expander.N = length(p_exp);
cycle_data.expander.eta = eta_expander;

% Pump thermodynamic trajectory
cycle_data.pump_h.p = p_pump_h;
cycle_data.pump_h.h = h_pump_h;
cycle_data.pump_h.T = T_pump_h;
cycle_data.pump_h.s = s_pump_h;
cycle_data.pump_h.d = d_pump_h;
cycle_data.pump_h.N = length(p_pump_h);
cycle_data.pump_h.eta = eta_pump_h;

% Pump thermodynamic trajectory
cycle_data.pump_f.p = p_pump_f;
cycle_data.pump_f.h = h_pump_f;
cycle_data.pump_f.T = T_pump_f;
cycle_data.pump_f.s = s_pump_f;
cycle_data.pump_f.d = d_pump_f;
cycle_data.pump_f.N = length(p_pump_f);
cycle_data.pump_f.eta = eta_pump_f;

% Pump thermodynamic trajectory
cycle_data.pump_c.p = p_pump_c;
cycle_data.pump_c.h = h_pump_c;
cycle_data.pump_c.T = T_pump_c;
cycle_data.pump_c.s = s_pump_c;
cycle_data.pump_c.d = d_pump_c;
cycle_data.pump_c.N = length(p_pump_c);
cycle_data.pump_c.eta = eta_pump_c;

% Evaporator thermodynamic trajectory
cycle_data.evaporator.p_h = p_evap_h;
cycle_data.evaporator.p_c = p_evap_c;
cycle_data.evaporator.T_h = T_evap_h;
cycle_data.evaporator.T_c = T_evap_c;
cycle_data.evaporator.h_h = h_evap_h;
cycle_data.evaporator.h_c = h_evap_c;
cycle_data.evaporator.d_h = d_evap_h;
cycle_data.evaporator.d_c = d_evap_c;
cycle_data.evaporator.s_h = s_evap_h;
cycle_data.evaporator.s_c = s_evap_c;
cycle_data.evaporator.dT = dT_evap;
cycle_data.evaporator.N = N_evap;

% Condenser thermodynamic trajectory
cycle_data.condenser.p_h = p_cond_h;
cycle_data.condenser.p_c = p_cond_c;
cycle_data.condenser.T_h = T_cond_h;
cycle_data.condenser.T_c = T_cond_c;
cycle_data.condenser.h_h = h_cond_h;
cycle_data.condenser.h_c = h_cond_c;
cycle_data.condenser.d_h = d_cond_h;
cycle_data.condenser.d_c = d_cond_c;
cycle_data.condenser.s_h = s_cond_h;
cycle_data.condenser.s_c = s_cond_c;
cycle_data.condenser.dT = dT_cond;
cycle_data.condenser.N = N_cond;

% Recuperator thermodynamic trajectory
cycle_data.recuperator.p_h = p_rec_h;
cycle_data.recuperator.p_c = p_rec_c;
cycle_data.recuperator.T_h = T_rec_h;
cycle_data.recuperator.T_c = T_rec_c;
cycle_data.recuperator.h_h = h_rec_h;
cycle_data.recuperator.h_c = h_rec_c;
cycle_data.recuperator.d_h = d_rec_h;
cycle_data.recuperator.d_c = d_rec_c;
cycle_data.recuperator.s_h = s_rec_h;
cycle_data.recuperator.s_c = s_rec_c;
cycle_data.recuperator.dT = dT_rec;
cycle_data.recuperator.N = N_rec;

% Fluid properties
cycle_data.properties.p_0 = p_0;
cycle_data.properties.T_0 = T_0;
cycle_data.properties.p_crit = p_crit;
cycle_data.properties.T_crit = T_crit;
cycle_data.properties.p_trip = p_trip;
cycle_data.properties.T_trip = T_trip;
cycle_data.properties.p_max = p_max;
cycle_data.properties.T_max = T_max;
cycle_data.properties.h_max = h_max;
cycle_data.properties.p_min = p_min;
cycle_data.properties.T_min = T_min;
cycle_data.properties.h_min = h_min;

% Energy and exergy check
cycle_data.checks.energy1 = energy_error1;
cycle_data.checks.energy2 = energy_error2;
cycle_data.checks.exergy = exergy_error;


%% Export the optimization problem variables
% Degrees of freedom
cycle_data.optimization.x = x;

% Lower bounds
cycle_data.optimization.lb = [fixed_parameters.lower_bounds.x1,   ...
                              fixed_parameters.lower_bounds.x2,   ...
                              fixed_parameters.lower_bounds.x3,   ...
                              fixed_parameters.lower_bounds.x4,   ...
                              fixed_parameters.lower_bounds.x5,   ...
                              fixed_parameters.lower_bounds.x6,   ...
                              fixed_parameters.lower_bounds.x7]';

% Upper bounds
cycle_data.optimization.ub = [fixed_parameters.upper_bounds.x1,   ...
                              fixed_parameters.upper_bounds.x2,   ...
                              fixed_parameters.upper_bounds.x3,   ...
                              fixed_parameters.upper_bounds.x4,   ...
                              fixed_parameters.upper_bounds.x5,   ...
                              fixed_parameters.upper_bounds.x6,   ...
                              fixed_parameters.upper_bounds.x7]';
                          
% Optimization constraints
cycle_data.optimization.c = c_ineq;
cycle_data.optimization.c_eq = c_eq;

% Store bounds information in a cell
bounds_cell{1,1} = 'Variable name';
bounds_cell{1,2} = 'Lower bound';
bounds_cell{1,3} = 'Variable value';
bounds_cell{1,4} = 'Upper bound'; k = 2;
bounds_cell(k,1) = {'Heat source exit temperature'}; k = k+1;
bounds_cell(k,1) = {'Heat sink exit temperature'}; k = k+1;
bounds_cell(k,1) = {'Pump inlet pressure'}; k = k+1;
bounds_cell(k,1) = {'Pump inlet enthalpy'}; k = k+1;
bounds_cell(k,1) = {'Expander inlet pressure'}; k = k+1;
bounds_cell(k,1) = {'Expander inlet enthalpy'}; k = k+1;
bounds_cell(k,1) = {'Recuperator effectiveness'};
bounds_cell(2:k,2) = num2cell(cycle_data.optimization.lb);
bounds_cell(2:k,3) = num2cell(cycle_data.optimization.x);
bounds_cell(2:k,4) = num2cell(cycle_data.optimization.ub);
cycle_data.optimization.check_bounds = bounds_cell;

% Store constraints information in a cell
% Headers
constraints_cell{1,1} = 'Constraint name';
constraints_cell{1,2} = 'Lower limit';
constraints_cell{1,3} = 'Constraint value';
constraints_cell{1,4} = 'Upper limit';
constraints_cell{1,5} = 'Applied?'; k = 2;

% Evaporator minimum temperature difference
constraints_cell(k,1) = {'Evaporator pinch point'};
constraints_cell(k,2) = {dT_evap_min};
constraints_cell(k,3) = {min(dT_evap)};
constraints_cell(k,4) = {inf};
constraints_cell(k,5) = {fixed_parameters.constraints.dT_evap.apply}; k = k+1;

% Condenser minimum temperature difference
constraints_cell(k,1) = {'Condenser pinch point'};
constraints_cell(k,2) = {dT_cond_min};
constraints_cell(k,3) = {min(dT_cond)};
constraints_cell(k,4) = {inf};
constraints_cell(k,5) = {fixed_parameters.constraints.dT_cond.apply}; k = k+1;

% Recuperator minimum temperature difference
constraints_cell(k,1) = {'Recuperator pinch point'};
constraints_cell(k,2) = {dT_rec_min};
constraints_cell(k,3) = {min(dT_rec)};
constraints_cell(k,4) = {inf};
constraints_cell(k,5) = {fixed_parameters.constraints.dT_rec.apply}; k = k+1;

% Subcooling constraint
constraints_cell(k,1) = {'Expander vapor quality'};
constraints_cell(k,2) = {expander_quality_min};
constraints_cell(k,3) = {min(q_exp)};
constraints_cell(k,4) = {inf};
constraints_cell(k,5) = {fixed_parameters.constraints.quality.apply}; k = k+1;

% Subcooling constraint
constraints_cell(k,1) = {'Degree of subcooling'};
constraints_cell(k,2) = {dT_subcooling_min};
constraints_cell(k,3) = {dT_subcooling};
constraints_cell(k,4) = {dT_subcooling_max};
constraints_cell(k,5) = {fixed_parameters.constraints.dT_subcooling.apply}; k = k+1;

% Superheating constraint
constraints_cell(k,1) = {'Degree of superheating'};
constraints_cell(k,2) = {dT_superheating_min};
constraints_cell(k,3) = {dT_superheating};
constraints_cell(k,4) = {dT_superheating_max};
constraints_cell(k,5) = {fixed_parameters.constraints.dT_superheating.apply}; k = k+1;

% Condenser minimum temperature difference
constraints_cell(k,1) = {'Cycle minimum pressure (bar)'};
constraints_cell(k,2) = {p_min_constraint/1e5};
constraints_cell(k,3) = {min(p_cycle)/1e5};
constraints_cell(k,4) = {p_max_constraint/1e5};
constraints_cell(k,5) = {fixed_parameters.constraints.pressure.apply}; k = k+1;

% Condenser minimum temperature difference
constraints_cell(k,1) = {'Cycle maximum pressure (bar)'};
constraints_cell(k,2) = {p_min_constraint/1e5};
constraints_cell(k,3) = {max(p_cycle)/1e5};
constraints_cell(k,4) = {p_max_constraint/1e5};
constraints_cell(k,5) = {fixed_parameters.constraints.pressure.apply};

% Store the cell with constraint check
cycle_data.optimization.check_constraints = constraints_cell;


end