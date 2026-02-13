% Spectral Synthesis of 2D random field with given correlation structure

function [Y,radial_spectra,DD] = SpecSyn2D_Hybrid(Deterministic, wavelength_cut,N, samp, corr, acf, Rseed)

% By Victor Hernández (victorh@hi.is). July 2025. 
% Based on code by Martin Mai's code SpecSyn2D

% Generates a 2D-random field Y of size (Nz+1 x Nx+1) with possibly variable spatial 
% sampling in z,x-direction. At wavenumbers < wavelength_cut the Deterministic slip is used.
% This function simulates anisotropic random fields, i.e.   
% the correlation len in both directions can be different rectangular grid dimensions 
% are also possible corr contains the corr. len ax, az and the Hurst number H or the fractal 
% dimension  D. The autocorrelation function to  be used has to be specified by acf
%  NOTE: D = 3-H, and larger D (smaller H) yield "rougher" fields
% The algorithm is based on the spectral synthesis method (Ogorodnikov and Prigarin, 1996)
% The merging of low and high wavenumber parts is according to Graves and Pitarka (2010)
% Ogorodnikov, V.A. and Prigarin, S.M. (1996), Numerical Modelling of Random Processes 
% and Fields: Algorithms and Applications. De Gruyter.
%  INPUT:
% Deterministic - Deterministic 2D slip spatial field covering the same  area as N
%  N     - grid dimensions [Nz Nx] in m
% samp   - desired sampling, [dz dx] if scalar, then dz = dx
% corr  - corr = [az ax]   for 'gs' and 'ex' 
%     corr = [az ax H] for 'ak' note that 0 <= H <= 1
%     corr = [D kc]    for 'fr' D is the fractal dimension,
%      kc: corner wavenumber, spectrum decays linearly for k>kc
%  acf   - autocorrelation function: 
%    'gs' or 'GS' - Gaussian
%    'ex' or 'EX' - Exponential
%    'ak' or 'AK' - anisotropic vonKarman
%    'fr' or 'FR' - fractal distribution
%  Rseed - seeds for random number generators if omitted or empty
%     Rseed = sum(100*clock) is used (returned in structure spar)
%     [Rpseed Rsseed] for the phase and small random spectral part

% OUTPUT: 2D-random field of size (Nz+1 x Nx+1) approximately standard Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6, acf = 'ak'; end
if nargin < 7, Rseed = []; end

if length(samp) == 1, samp = [samp samp]; end
if length(N) == 1, N = [N N]; end

if mod(N(2), samp(2)) ~= 0
    error('Sampling in X does not yield integer grid points. Abort.');
end
if mod(N(1), samp(1)) ~= 0
    error('Sampling in Z does not yield integer grid points. Abort.');
end

if strcmpi(acf,'fr')
    if length(corr) == 2
        D = corr(1); kc = corr(2);
    else
        D = corr(1); kc = 0.1;
    end
elseif strcmpi(acf,'ex') || strcmpi(acf,'gs')
    az = corr(1); ax = corr(2);
elseif strcmpi(acf,'ak')
    az = corr(1); ax = corr(2); H = corr(3); 
end

% set size for spectral synthesis tool that generates fields of size (2*rmz+1, 2*rmx+1), 
% i.e. the method requires an ODD number of points in each direction
nptsX = N(2)/samp(2);  %% number of grid-points in X
nptsZ = N(1)/samp(1);  %% number of grid-points in Z

if     mod(nptsX,2) == 0, rmx = nptsX/2;
elseif mod(nptsX,2) == 1, rmx = (nptsX-1)/2;
end
if     mod(nptsZ,2) == 0, rmz = nptsZ/2;
elseif mod(nptsZ,2) == 1, rmz = (nptsZ-1)/2;
end

kx = (rmx./((2*rmx+1)*samp(2)))*linspace(-2*pi,2*pi,2*rmx+1);
kz = (rmz./((2*rmz+1)*samp(1)))*linspace(-2*pi,2*pi,2*rmz+1);
%%% compose power spectrum for two of the four quadrants
kr = zeros(rmz+1,2*rmx+1);

for j = 1:(2*rmx+1)
    for i = 1:(rmz+1)
	if ((acf == 'fr') | (acf == 'FR'))
	   kr(i,j) = sqrt((kz(i)^2) + (kx(j)^2));
	else
	   kr(i,j) = sqrt((az^2)*(kz(i)^2) + (ax^2)*(kx(j)^2));
	end
    end
end

% Power Spectrum
switch lower(acf)
    case 'gs'
        PS = 0.5 * ax * az * exp(-0.25 * kr.^2);
    case 'ex'
        PS = (ax * az) ./ (1 + kr.^2).^1.5;
    case 'ak'
        % coef = 4*pi*H*ax*az / besselk(H,min(kr(kr>0)));
        coef = 1; %does not matter
        PS = coef ./ (1 + kr.^2).^(H+1);
    case 'fr'
        decay = 0.5*(8 - 2*D);
        kr(kr==0) = pi*kc;
        PS = 1 ./ (kr.^2).^decay;
end

if isempty(Rseed)
    rng('shuffle');
else
    rng(Rseed(1));
end

% Deterministic slip for low wavenumbers DD<->DD_fft
[Nz_D,Nx_D]=size(Deterministic);
[x_D,z_D]=meshgrid(linspace(0,N(2),Nx_D),linspace(0,N(1),Nz_D));
[x,z]=meshgrid(linspace(0,N(2),2*rmx+1),linspace(0,N(1),2*rmz+1));
DD=interp2(x_D,z_D,Deterministic,x,z,'linear');
DD_fft = fft2(DD);

% Filter according to Graves and Pitarka (2010)
[KX, KZ] = meshgrid(ifftshift(kx),ifftshift(kz));
kf = sqrt((wavelength_cut/2/pi)^2*KX.^2 + (wavelength_cut/2/pi)^2*KZ.^2);
kr2 = sqrt(KX.^2 + KZ.^2);
filter=1./(1 + kf.^8);

% Generating stochastic part
[m,n] = size(kr);
xi1 = randn(m,n); xi2 = randn(m,n);
SpectralSynthesis = sqrt(0.5*PS).* (xi1 + 1i*xi2); %cartesian form
SpectralSynthesis(rmz+1,rmx+1) = sqrt(PS(rmz+1,rmx+1)) * 1; %*randn
D_SpectralSynthesis=HermitianSymmetry(SpectralSynthesis,rmx,rmz);

% Scaling factor for the stochastic part
num = sum(abs(DD_fft(kr2<2*pi/wavelength_cut)), 'all'); 
den = sum(abs(D_SpectralSynthesis(kr2<2*pi/wavelength_cut)), 'all');
s   = max(num,0) / max(den, eps);
D_SpectralSynthesis=D_SpectralSynthesis.*s;

% Merging
spectrum = DD_fft.*filter +D_SpectralSynthesis .*(1-filter);

% Inverse FFT
Y = real(ifft2(spectrum));
Y = Y / std(Y(:)); % to have unit variance

if mean(Y(:)) < 0; Y = (-1)*Y; end	% positive mean (though small)

%Resampling to get required number of points
[x_out,z_out]=meshgrid(linspace(0,N(2),nptsX),linspace(0,N(1),nptsZ));
Y = interp2(x,z,Y,x_out,z_out,'linear');
DD = interp2(x,z,DD,x_out,z_out,'linear');

%% Spectra check for debugging - Not necessary
%  %[krad_f, Pkrad_f] = radial_avg_2D(kr2, 1-filter); figure; loglog(krad_f,Pkrad_f); hold on; plot([2*pi/wavelength_cut 2*pi/wavelength_cut],[1E0 1e-2],'r--')
 %Stochastic spectrum;
 [krad, Pkrad_S] = radial_avg_2D(kr2, abs(D_SpectralSynthesis));
 %checking the decay of the spectrum
 % P = polyfit(log10(krad(35:45)),log10(Pkrad_S(35:45)),1);P(1); 
 %Deterministic spectrum
 [krad_D, Pkrad_D] = radial_avg_2D(kr2, abs(DD_fft));% figure; loglog(krad_D,Pkrad_D)
 % P = polyfit(log10(krad(18:30)),log10(Pkrad_D(18:30)),1);P(1);
%Hybrid
 [krad_T, Pkrad_T] = radial_avg_2D(kr2, abs(spectrum)); 
%  figure; loglog(krad,Pkrad_S,'k-',krad_D,Pkrad_D,'b-',krad_T,Pkrad_T,'g--','linewidth',1.5);hold on; plot([2*pi/wavelength_cut 2*pi/wavelength_cut],[max(Pkrad_T) max(Pkrad_T)/100],'r--')
%  legend('Stochastic','Deterministic','Merged','k_c','FontSize',12)
%  grid on
%  set(gca,'FontSize',12)
%  ylim([1,1e4])
%  xlabel('Radial wavenumber [rad/m]'); ylabel('Spectral amplitude')
radial_spectra.k=krad; radial_spectra.Pkrad_S=Pkrad_S;radial_spectra.Pkrad_D=Pkrad_D;radial_spectra.Pkrad_T=Pkrad_T;

end

%%
function [U] = HermitianSymmetry(Y,rmx,rmz)
U = zeros(2*rmz+1,2*rmx+1);	        %% will be conj. sym. field
Y = [Y ; conj(fliplr(flipud(Y(1:rmz,:))))];

for i = 1:rmx
     Y(rmz+1,-i+2*rmx+2) = conj(Y(rmz+1,i));
end
for i = 1:rmz+1
   for j = 1:rmx+1
      U(i,j) = Y(i+rmz,j+rmx);
    end
end
for i = rmz+2:2*rmz+1
   for j = rmx+2:2*rmx+1
      U(i,j) = Y(i-rmz-1,j-rmx-1);
    end
end
for i = 1:rmz+1
   for j = rmx+2:2*rmx+1
      U(i,j) = Y(i+rmz,j-rmx-1);
    end
end
for i = rmz+2:2*rmz+1
   for j = 1:rmx+1
      U(i,j) = Y(i-1-rmz,j+rmx);
    end
end

end
