function q = quality(prop2,prop2_value,fluid,p,p_crit)

% Give an output for quality even if the pressure is supercritical
if p < p_crit
    q = prop_calculation('Q','P',p,prop2,prop2_value,fluid);
elseif p>= p_crit
    q = 1.1;
else
    error('Ooops, something went wrong in quality computation');
end

if q < 0
    q = 1.1;

end