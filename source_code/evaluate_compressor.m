function [h_high,p,h,T,s,d] = evaluate_compressor(fluid,h_low,p_low,p_high,eta_poly,calc_detail)

% Compute the thermodynamic trajectory
% options = odeset('RelTol',1e-6,'AbsTol',1e-6);
% [p,h] = ode45(@(p,h)R(p,h,fluid,eta_poly),[p_low p_high],h_low,options);
% h_high = h(end);

s = prop_calculation('S','P',p_low,'H',h_low,fluid);
h_high_s = prop_calculation('H','P',p_high,'S',s,fluid);
h_high = h_low+(h_high_s-h_low)/eta_poly;
p = [p_low p_high];
h = [h_low h_high];

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

% function dhdp = R(p,h,fluid,eta_poly)
% rho = prop_calculation('d','p',p,'h',h,fluid);
% dhdp = 1000/rho/eta_poly;
% end


% N = 10;
% dp = (p_high-p_low)/(N-1);
% p = (p_low:dp:p_high)';
% h = zeros(N,1);
% h(1) = h_low;
% for n = 1:N-1
%     h(n+1) = h(n)+dp*R(p(n),h(n),fluid,eta_poly);
% end
% h_high = h(end);