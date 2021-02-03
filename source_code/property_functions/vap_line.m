function [X_vap, Y_vap] = vap_line(fluid,T_min,T_max,x,y,N_vap)

% Clustering parameter (>1)
beta = 1.0010;               

% Vapor saturation line
z_vap = linspace(0,1,N_vap)';
T_vap = T_max+cluster_func(z_vap,beta)*(T_min-T_max);
X_vap = zeros(N_vap,1);                 % Pre-allocate space
Y_vap = zeros(N_vap,1);                 % Pre-allocate space
for i = 1:N_vap    
    X_vap(i) = prop_calculation(x,'T',T_vap(i),'Q',1,fluid);    
    Y_vap(i) = prop_calculation(y,'T',T_vap(i),'Q',1,fluid);
end

end