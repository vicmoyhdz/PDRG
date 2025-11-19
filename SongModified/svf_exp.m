% Coded by Song, S. (April 2010)
% Generating the shape of slip velocity function (SVF)
% based on Liu et al. (BSSA, 2006)

function [svf, t, sf] = svf_exp(tau,dt,nt)
tau=tau/6.64;
t = [0:dt:dt*(nt-1)];
%%% End of inputs

svf = (t/tau^2).*exp(-t/tau);
sf =1-(t/tau+1).*exp(-t/tau);

end
