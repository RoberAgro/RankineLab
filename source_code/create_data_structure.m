function cycle_data = create_data_structure

% This function is not really necessary for the computations since it is
% not necessary to declare variables in MATLAB
% It can be used as a reference where the cycle variables are documented

% States of the thermodynamic cycle
cycle_states = struct('fluid',       [], ...                               % Fluid name
                      'p',           [], ...                               % Pressure [kPa]
                      'T',           [], ...                               % Temperature [K]
                      'd',           [], ...                               % Density [kg/m3]
                      'h',           [], ...                               % Enthalpy [J/kg]
                      's',           [], ...                               % Entropy [J/kg K]
                      'e',           []);                                  % Exergy [J/kg]
                  
cycle_states(1:7) = cycle_states;                                          

%Fluid properties
properties   = struct('p_0',           [], ...                             % Ambient pressure
                      'T_0',           [], ...                             % Ambient temperature
                      'p_crit',        [], ...                             % Critical pressure
                      'T_crit',        [], ...                             % Critical temperature
                      'p_trip',        [], ...                             % Triple pressure
                      'T_trip',        []);                                % Triple temperature

% Working fluids
fluids = [];
     

% Mass flows
mass_flows   = struct('m_h',         [], ...                               % Heating fluid mass flow rate [kg/s]
                      'm_f',         [], ...                               % Working fluid mass flow rate [kg/s]
                      'm_c',         []);                                  % Cooling fluid mass flow rate [kg/s]

                  
% Energy flows
energy_flows = struct('Q_max',        [], ...                              % Heat flow rate if the heat source is cooled down to ambient temperature [W]
                      'Q_evap',       [], ...                              % Heat flow rate absorbed in the heat exchanger [W]
                      'W_exp',        [], ...                              % Power delivered by the expander [W]
                      'W_pump_h',     [], ...                              % Power consumed by the pump [W]
                      'W_pump_f',     [], ...                              % Power consumed by the pump [W]
                      'W_pump_c',     []);                                 % Power consumed by the pump [W]
          

% % Exergy flows
exergy_flows = struct('E_max',       [], ...                               % Exergy flow if the heat source is cooled down to ambient temperature [W]
                      'I_evap',      [], ...                               % Exergy destruction in the heat exchanger [W]
                      'I_cond',      [], ...                               % Exergy destruction in the heat exchanger [W]
                      'I_rec',       [], ...                               % Exergy destruction in the heat exchanger [W]
                      'I_pump_h',    [], ...                               % Exergy destruction in the pump [W]
                      'I_pump_f',    [], ...                               % Exergy destruction in the pump [W]
                      'I_pump_c',    [], ...                               % Exergy destruction in the pump [W]
                      'I_exp',       [], ...                               % Exergy destruction in the expander [W]
                      'I_source',    [], ...                               % Exergy flow at the exit of the heat source [W]
                      'I_sink',      []);                                  % Exergy flow at the exit of the heat sink [W]   


% Efficiencies of the cycle
efficiency = struct('eta_1',         [], ...                               % Plant energy efficiency
                    'eta_1recovery', [], ...                               % Energy recovery efficiency
                    'eta_1cycle',    [], ...                               % Cycle energy efficiency
                    'eta_2',         [], ...                               % Plant exergy efficiency
                    'eta_2recovery', [], ...                               % Exergy recovery efficiency
                    'eta_2cycle',    []);                                  % Cycle exergy efficiency

                  
% Variables associated with the expander
expander     = struct('N',            [], ...                              % Number of discretizations
                      'p',            [], ...                              % Pressure [kPa]
                      'T',            [], ...                              % Temperature [K]
                      'd',            [], ...                              % Density [kg/m3]
                      'q',            [], ...                              % Quality [-]
                      'h',            [], ...                              % Enthalpy [J/kg]
                      's',            [], ...                              % Entropy [J/kg K]
                      'eta',          []);                                 % Polytropic efficiency
                  
% Variables associated with the pump
pump         = struct('N',           [], ...                               % Number of discretizations
                      'p',            [], ...                              % Pressure [kPa]
                      'T',            [], ...                              % Temperature [K]
                      'd',            [], ...                              % Density [kg/m3]
                      'q',            [], ...                              % Quality [-]
                      'h',            [], ...                              % Enthalpy [J/kg]
                      's',            [], ...                              % Entropy [J/kg K]
                      'eta',          []);                                 % Polytropic efficiency
                  
          
% Variables associated with the heat exchanger
exchanger   = struct('N',            [], ...                               % Number of discretizations
                     'p_h',          [], ...                               % Pressure [kPa]
                     'p_c',          [], ...                               % Pressure [kPa]
                     'T_h',          [], ...                               % Temperature [K]
                     'T_c',          [], ...                               % Temperature [K]
                     'd_h',          [], ...                               % Density [kg/m3]
                     'd_c',          [], ...                               % Density [kg/m3]
                     'h_h',          [], ...                               % Enthalpy [J/kg]
                     'h_c',          [], ...                               % Enthalpy [J/kg]
                     's_h',          [], ...                               % Entropy [J/kg K]
                     's_c',          [], ...                               % Entropy [J/kg K]
                     'dT_min',       []);                                  % Pinch point [K]
                  

% Variables associated with the optimization problem
optimization = struct('x',                   [], ...                       % Vector of degrees of freedom
                      'lb',                  [], ...                       % Vector of upper bounds
                      'ub',                  [], ...                       % Vector of lower bounds
                      'c',                   [], ...                       % Vector of inequality constraints
                      'c_eq',                []);                          % Vector of equality constraints
 


% Create the cycle_data structure grouping the sub-structures
cycle_data = struct('cycle_states',       cycle_states, ...
                    'properties',         properties,  ...
                    'fluids',             fluids,        ...
                    'mass_flows',         mass_flows,    ...
                    'energy_flows',       energy_flows,  ...
                    'exergy_flows',       exergy_flows,  ...     
                    'efficiency',         efficiency,    ...
                    'evaporator',         exchanger,     ...
                    'condenser',          exchanger,     ...
                    'recuperator',        exchanger,     ...
                    'expander',           expander,       ...
                    'pump_h',             pump,          ...
                    'pump_f',             pump,          ...
                    'pump_c',             pump,          ...
                    'optimization',       optimization);
                      

end
