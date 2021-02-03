function [dT,T_h,T_c,h_h,h_c,p_h,p_c,d_h,d_c,s_h,s_c] = evaluate_exchanger(fluid_h,fluid_c,m_h,m_c,p_h1,p_h2,p_c1,p_c2,h_h1,h_h2,h_c1,h_c2,N,calc_detail)
%% Preliminary computations

% Compute the mass flow rate ratio
m = m_h/m_c;

% Compute the temperature difference at the inlet
T_h1 = prop_calculation('T','P',p_h1,'H',h_h1,fluid_h);
T_c1 = prop_calculation('T','P',p_c1,'H',h_c1,fluid_c);


%% Discretize the heat exchanger in N nodes
% Discretize the heat exchanger in N nodes
dq = (h_h2-h_h1)/(N-1);     % Heat exchange step
dp_h = (p_h2-p_h1)/(N-1);   % Hot pressure drop step
dp_c = (p_c2-p_c1)/(N-1);   % Cold pressure drop step

% Preallocate memory
h_h = zeros(N,1); h_c = zeros(N,1);
T_h = zeros(N,1); T_c = zeros(N,1);
p_h = zeros(N,1); p_c = zeros(N,1);

% Initialize variables
h_h(1) = h_h1; h_c(1) = h_c1;
T_h(1) = T_h1; T_c(1) = T_c1;
p_h(1) = p_h1; p_c(1) = p_c1;

% Main loop
for i = 1:N-1
    h_h(i+1) = h_h(i)+dq;
    h_c(i+1) = h_c(i)+dq*m;
    p_h(i+1) = p_h(i)+dp_h;
    p_c(i+1) = p_c(i)+dp_c;
    T_h(i+1) = prop_calculation('T','P',p_h(i+1),'H',h_h(i+1),fluid_h);
    T_c(i+1) = prop_calculation('T','P',p_c(i+1),'H',h_c(i+1),fluid_c);
%     q_h(i+1) = prop_calculation('Q','P',p_h(i+1),'H',h_h(i+1),fluid_h);
%     q_c(i+1) = prop_calculation('Q','P',p_c(i+1),'H',h_c(i+1),fluid_c);
end

% I am not satisfied with the solution adding extra nodes at phase changes
% The linear pressure drop with heat transfer is kept to have a close
% thermodynamic diagram when plotting
% I think that the extra programming complexity and loss of robustness is
% not worth the extra accuracy at the phase changes
% If an accurate cycle simulation is desired the number of the
% discretizations can be increased after a local optima has been found


%% Compute the other thermodynamic variables
if strcmp(calc_detail,'short')
    d_h = zeros(N,1); d_c = zeros(N,1); s_h = zeros(N,1); s_c = zeros(N,1); 
elseif strcmp(calc_detail,'long')
    d_h = zeros(N,1); d_c = zeros(N,1); s_h = zeros(N,1); s_c = zeros(N,1); 
    for i = 1:N
        d_h(i) = prop_calculation('D','P',p_h(i),'H',h_h(i),fluid_h);
        d_c(i) = prop_calculation('D','P',p_c(i),'H',h_c(i),fluid_c);
        s_h(i) = prop_calculation('S','P',p_h(i),'H',h_h(i),fluid_h);        
        s_c(i) = prop_calculation('S','P',p_c(i),'H',h_c(i),fluid_c);
    end
end


%% Compute the pinch point temperature difference
dT = T_h-T_c;
% dT_min = min(dT);


end
