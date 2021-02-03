function [X_sat, Y_sat] = sat_line(fluid,T_min,T_max,x,y,N)

% Compute the liquid and the vapor saturation lines
[X_liq, Y_liq] = liq_line(fluid,T_min,T_max,x,y,ceil(N/2));
[X_vap, Y_vap] = vap_line(fluid,T_min,T_max,x,y,floor(N/2));

% Export the results
X_sat = [X_liq; X_vap];    % Create a single vector with property x
Y_sat = [Y_liq; Y_vap];    % Create a single vector with property y

end