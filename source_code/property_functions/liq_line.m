function [X_liq, Y_liq] = liq_line(fluid,T_min,T_max,x,y,N_liq)

% Clustering parameter (>1)
beta = 1.0010;

% Liquid saturation line
z_liq = linspace(0,1,N_liq)';
T_liq = T_max+cluster_func(z_liq,beta)*(T_min-T_max);
T_liq(1:end) = T_liq(end:-1:1);         % Reverse the vector (clustering)
X_liq = zeros(N_liq,1);                 % Pre-allocate space
Y_liq = zeros(N_liq,1);                 % Pre-allocate space
for i = 1:N_liq
    X_liq(i) = prop_calculation(x,'T',T_liq(i),'Q',0,fluid);
    Y_liq(i) = prop_calculation(y,'T',T_liq(i),'Q',0,fluid);    
end

end