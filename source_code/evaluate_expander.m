function [h_low,p,h,T,s,d] = evaluate_expander(fluid,h_high,p_high,p_low,eta_expander,eta_expander_definition,calc_detail)

if strcmp(eta_expander_definition,'isentropic') == 1
    
    % Isentropic efficiency computation
    s = prop_calculation('S','P',p_high,'H',h_high,fluid);
    h_low_s = prop_calculation('H','P',p_low,'S',s,fluid);
    h_low = h_high-(h_high-h_low_s)*eta_expander;
    p = [p_low p_high]';
    h = [h_low h_high]';
    
elseif strcmp(eta_expander_definition,'polytropic') == 1
    
    % Polytropic efficiency computation
    options = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [p,h] = ode45(@(p,h)R(p,h,fluid,eta_expander),linspace(p_high,p_low,25),h_high,options);
    h_low = h(end);
    
else
    
    error("The expander efficiency definition must be 'isentropic' or 'polytropic'")
    
end
    
% Compute the other thermodynamic variables
if strcmp(calc_detail,'short')
    T = 0*p; s = 0*p; d = 0*p;
elseif strcmp(calc_detail,'long')
    T = 0*p; s = 0*p; d = 0*p;
    for i = 1:length(p)
        T(i) = prop_calculation('T','P',p(i),'H',h(i),fluid);
        s(i) = prop_calculation('S','P',p(i),'H',h(i),fluid);
        d(i) = prop_calculation('D','P',p(i),'H',h(i),fluid);
    end
end

end

function dhdp = R(p,h,fluid,eta_poly)
rho = prop_calculation('D','P',p,'H',h,fluid);
dhdp = eta_poly/rho;
end

% N = 10;
% dp = (p_low-p_high)/(N-1);
% p = (p_high:dp:p_low)';
% h = zeros(N,1);
% h(1) = h_high;
% for n = 1:N-1
%     h(n+1) = h(n)+dp*R(p(n),h(n),fluid,eta_poly);
% end
% h_low = h(end);