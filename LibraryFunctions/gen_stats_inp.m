%%% Coded by Song, S.G. (June 2015)
%close all, clear all

function [rup] = gen_stats_inp(rup)

ifl = 'src_stats_mod.mat'; % reading source statistics model

A = rup.L*rup.W;
mu_slip = rup.target_Mo/rup.mu/A; % meter to cm

load(ifl) 
mu_mu(3)=rup.user_stats.mean_Vmax;
mu_mu(2)=rup.user_stats.mean_Vr;

mu_sig(3)=rup.user_stats.sigma_Vmax;
mu_sig(2)=rup.user_stats.sigma_Vr;

L = chol(C,'lower');

N = length(C);

rng(rup.seed.seed,'twister'); 
for i=1:rup.Nstats
  
  s(:,i) = randn(N,1);
  
  s(1,i) = (mu_slip-mu_mu(1))/sig_mu(1)/L(1);  %Conditioning on a specific value of slip

  s(:,i) = L*s(:,i);

  s(:,i) = s(:,i).*[sig_mu sig_sig sig_ax sig_az sig_cc]' ...
      + [mu_mu mu_sig mu_ax mu_az mu_cc]';

  s(7:18,i) = exp(s(7:18,i));  %Transforming to Lognormal distribution
  
end

src_stats_inp = s'; 

rup.stats = src_stats_inp;



