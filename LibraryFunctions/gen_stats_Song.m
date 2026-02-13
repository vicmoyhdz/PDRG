%%% Coded by Song, S. (April. 2013)
%%% Set up target 1-point and 2-point statistics

function [rup] = gen_stats_Song(rup)

rup.Nstats = 1000;
rup.stats_id = 100;

%%% 1-point and 2-point statistics of source parameters

%%%Sampling input 1-point and 2-point statistics, given the source statistics model
[rup] = gen_stats_inp(rup);

rup.p1.mu =  [mean(rup.slip_nontaper(:)) rup.user_stats.mean_Vr rup.user_stats.mean_Vmax];
rup.p1.sig =  [std(rup.slip_nontaper(:)) rup.user_stats.sigma_Vr rup.user_stats.sigma_Vmax];


rup.p2.ax = [rup.stats(rup.stats_id,7) rup.stats(rup.stats_id,8)  rup.stats(rup.stats_id,9); ...
             nan                       rup.stats(rup.stats_id,10) rup.stats(rup.stats_id,11); ...
             nan                       nan                        rup.stats(rup.stats_id,12)];

rup.p2.az = [rup.stats(rup.stats_id,13) rup.stats(rup.stats_id,14) rup.stats(rup.stats_id,15); ...
             nan                        rup.stats(rup.stats_id,16) rup.stats(rup.stats_id,17); ...
             nan                        nan                        rup.stats(rup.stats_id,18)];


rup.p2.cc = [1    rup.stats(rup.stats_id,19) rup.stats(rup.stats_id,20); ...
             nan  1                          rup.stats(rup.stats_id,21); ...
             nan  nan                        1];

end
