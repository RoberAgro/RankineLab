function f = cluster_func(z,beta)
% Muller introduced this function in his notes (I should find a reference)
% 0<z<1 is an uniformly spaced vector
% beta>1 is the clustering parameter

f = 1 + beta*(1-((beta+1)/(beta-1)).^(1-z))./(1+((beta+1)/(beta-1)).^(1-z));

% The cluster function puts more points in the region where z is close to 0
% If you want more points in the region where z is close to 1 define the
% problem in the opposite manner and then reverse the resulting vector

end