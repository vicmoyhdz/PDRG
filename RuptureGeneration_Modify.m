% By Victor M. Hernández-Aguirre (victor@hi.is). July 2025
% A kinematic rupture generator for ground-motion simulations: 
% Validation and scenarios in South Iceland

clc
clear
close all
rng('shuffle')

addpath '.\LibraryFunctions'

%% Input and preprocessing

%General options
rup.Filename='TestDeterministic';
% A folder named Filename is created to save the model (rup.m) and plots
folder_name=[cd,'\',rup.Filename];
if ~isfolder(folder_name)
    mkdir (['./',rup.Filename]);
end
rup.target_Mw = 6.6; 
seed=100*rand; 
rup.sampling=[300 300]; %subfault size in m
Mo = 10^(1.5*rup.target_Mw + 9.05);
% Width and length of the fault. A scaling law could be used
% but for historical earthquake one can use deterministic dimensions (in m)
WL=[9600 39900]; %(in m)
rup.L = WL(2);  rup.W = WL(1); 
rup.stk  = 315; rup.dip  = 60; rup.rake  = 270; 
rup.mech='ds'; %ds: dip-slip %ss: for strike-slip
rup.rake_cov=0; %to add a slight variation to rake, std in degrees (0 no variation)
rup.ref = [530016.7664	4505390.882	-2000]; %upper right corner of the fault, needed to define geometry
%file containing the layering of the medium and their properties: vp,vs,p,thickness,qs,qp
rup.mech_file=[cd,'\','example_mech_prop.dat']; 
rup.z_ref=2100; % Top z coordinate of the mech file, i.e., ground surface z coordinate
alpha_roughness=0.000; %usually between 0.002-0.01. Set 0 to consider a planar fault
%--------------------------------------------------------------------------------
%Options for slip generation
rup.ax=10^(-2.95+0.61*rup.target_Mw)*1000; %correlation lengths (see Mai and Beroza, 2002)
rup.az=10^(-1.4+0.333*rup.target_Mw)*1000; 
rup.Hurst=0.75;
% ax=17.7 +0.35 *0.7*WL(2)/1000; % Melgar and Hayes (2019)
% az=6.7+0.41*0.7*WL(1)/1000;
cor_len=[rup.az, rup.ax, rup.Hurst]; %correlation lengths in m and Hurst exponent (Melgar 0.38)
rup.acf_type='ak'; %Von Karman: ak
rup.mu=2e10; %average shear modulus in Pa (to scale slip to magnitude)
rho=2500; Vs=sqrt(rup.mu /rho); Vp=1.8*Vs; lambda_ave=rho*(Vp^2-2*Vs^2);
taper_width=[2000,0,2000];     % Taper width for slip in m  [left/right, top, bottom]
taper_function = 'hn';
%---------------------------------------------------------------------------------
%Options in case a slip distributions is to be modified
wavelength_cut=500; % transition wavelength between stochastic and deterministic (in m)
%The file consists of a matrix with slip values in meters (rows along dip, colums along strike)
Deterministic_slip=load('Example_Deterministic_Slip.srcmod');
max_slip=4.0; % same as deterministic slip
%--------------------------------------------------------------------------------
%Options for kinematic rupture generation

% Savran and Olsen (2020). Kinematic Rupture Generator Based on 3‐D Spontaneous
%Rupture Simulations Along Geometrically Rough Faults
rup.method='Savran';  
% Song et al. (2013, 2014, 2016). Developing a generalized pseudo-dynamic source 
% model of Mw 6.5–7.0 to simulate strong ground motions
% rup.method='Song'; 
rup.num = 1;         % number of simulations
rup.user_stats.mean_Vmax=2; rup.user_stats.sigma_Vmax=1.2; 
rup.user_stats.mean_Vr=0.8; rup.user_stats.sigma_Vr=0.12; %Not used
 %Parameters of the Levy (stable) distribution for Vr. Can be tuned. Delta is the mode of the pdf
rup.user_stats.alpha=1.4;rup.user_stats.beta=-0.95;rup.user_stats.gamma=0.04;rup.user_stats.delta=0.84; 
%min/max values [slip(m) Vr/Vs(km/s) Vmax(m/s) risT(s) peak_time(s)];
rup.p1.min = [0   0.3  0.05  0.05 0.07];
rup.p1.max = [7  0.95   8      5   0.20];
%Liu et al. (2006): pliu; Tinti et al. (2005): etinti; Exponential: exp
rup.svf.type = 'pliu'; 
rup.dx1 =400; % grid size for sub-calculation (Chol, Eig)
rup.seed.seed = seed; rup.seed.stats = seed;
%Parameters for peak time distribution when using Liu: mean, std and
%correlation coefficient (rho_PT) with Vmax
rup.user_stats.mean_PT=0.11; rup.user_stats.sigma_PT=0.03; rup.user_stats.rho_PT=-0.4;

%--------------------------------------------------------------------------------
% Preprocessing velocity model and fault geometry
rup=preprocessing_fault(rup);

%% Slip distribution

% Generation of a slip distribution using the unconditional Fourier-based method for 
% stochastic fields simulation following a specified correlation structure.
% Following mostly the proposal of Aquib et al. (2025).
% In short, uses inverse Fourier transform of a random complex spectrum in wave-number 
% domain with a prescribed PSD (von Kármán).  

%----------------------------------------------------------------------------------------------
%Option 2. Add short wavelength noise to a given slip distribution
% 
[slip_nontaper, slip,slipp,radial_spectra,Deterministic_resampled] = CreateSlipModified(Deterministic_slip,wavelength_cut,max_slip,WL, rup.target_Mw, rup.mu_vec, rup.sampling, cor_len, rup.acf_type, taper_width, taper_function,seed);
rup.slip_deterministic=slip;
rup.slip_nontaper=slip_nontaper;
 %----------------------------------------------------------------
 %Compute seismic moments and disturbed rakes at each subfault
  rup.Mo_vec = rup.mu_vec.*rup.slip_deterministic(:).*(rup.sampling(2)*rup.sampling(1));
   [rake_srcmod, corr_coeff] = generate_correlated_field(rup.slip_deterministic,0.1, rup.rake_cov,rup.rake);  
  rup.rake_vec=rake_srcmod(:);

%% Hypocenter location

%Option 1. Probabilistic location (Mai et al., 2005)

% [slip_pos,shyp,dhyp,rup.shyp,rup.dhyp] = Stat_Hypo(rup);
% 
%-------------------------------------------------------------
%Option 2. Deterministic location

rup.shyp = 5000;      % along strike location (from top center) of hypocenter in m
rup.dhyp = 6928;        % along dip location (from top edge) of hypocenter in m
%-------------------------------------------------------------
 %Hypo coordinates in global coordinates
 HypLoc = [rup.shyp -rup.dhyp*cosd(rup.dip) -rup.dhyp*sind(rup.dip)];    
 [HypGlo(1),HypGlo(2)] = find_rotated_coords(HypLoc(1),HypLoc(2),-(rup.stk-90));
 HypGlo(1) = HypGlo(1) + rup.ref(1); HypGlo(2) = HypGlo(2) + rup.ref(2); HypGlo(3) = HypLoc(3) + rup.ref(3);
 rup.hypo_coor = HypGlo;

%% Vmax and Vrup from Covariance model

% Generates multivariate correlated Gaussian fields using Cholesky 
% factorization method, conditioned on given slip. Then modifies them 
% according to 1-point statistics. Afterwards computes Trup and Tau,
rup = gen_rup(rup); 

%% Calculation of roughness from integrated stress field

%Stress drop in strike and dip directions from Ripperger and Mai (2004)
[sigmaS,sigmaD] = Slip2Stress(slip,rup.rake,rup.sampling,rup.mu,lambda_ave,0.5);

%Integrated stress field from finding the scalar potential field whose
%gradient is the stress drop vector field. This solves the Poisson equation:
ISF2 = integrate_vector_field(sigmaS, sigmaD, rup.sampling);

ISF_norm=ISF2./max(ISF2,[],"all");
rup.ISF=ISF_norm;

% Compute roughness and normal and slip vectors
[Roughness,norm_vectors] = fault_roughness_frac(-ISF_norm,alpha_roughness,rup.L,rup.W,rup.sampling(2), rup.sampling(1));
rup.normal_vectors = rotate_vector_fault2global(norm_vectors,  rup.stk, rup.dip);
rup.slip_vectors = slip_vec_from_normal(rup.normal_vectors,rup.rake_vec);
rup.roughness=Roughness(:);

%Update coordinates considering roughness
rup.glo_coor(:,1)= rup.glo_coor(:,1)+rup.roughness.* rup.normal_vectors(:,1);
rup.glo_coor(:,2)= rup.glo_coor(:,2)+rup.roughness.* rup.normal_vectors(:,2);
rup.glo_coor(:,3)= rup.glo_coor(:,3)+rup.roughness.* rup.normal_vectors(:,3);


%% Export source MATE file in SPEED format 
% rup=export_MATE(rup);

%% Save
save([folder_name,'\','rup.mat'], 'rup')

%% Plotting
Plotting_source_parameters_merged
