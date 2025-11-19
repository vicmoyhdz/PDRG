function [Y] = SpecSyn2D(N, samp, corr, acf, Rseed)

% Slight modification from Martin Mai's code by Victor Hernández (vmh5@hi.is). July 2025

% Generates a 2D-random field Y of size (Nz+1 x Nx+1) with possibly variable spatial 
% sampling in z,x-direction. This function simulates anisotropic random fields, i.e.   
% the correlation len in both directions can be different rectangular grid dimensions 
% are also possible corr contains the corr. len ax, az and the Hurst number H or the fractal 
% dimension  D. The autocorrelation function to  be used has to be specified by acf
%  NOTE: D = 3-H, and larger D (smaller H) yield "rougher" fields
% The algorithm is based on the spectral synthesis method (Ogorodniko and Prigarin, 1996)
% Ogorodnikov, V.A. and Prigarin, S.M. (1996), Numerical Modelling of Random Processes 
% and Fields: Algorithms and Applications. De Gruyter.

% INPUT:
% N     - grid dimensions [Nz Nx] in m
% samp   - desired sampling, [dz dx] if scalar, then dz = dx
% corr  - corr = [az ax]   for 'gs' and 'ex' 
% corr = [az ax H] for 'ak' note that 0 <= H <= 1
% corr = [D kc]    for 'fr' D is the fractal dimension, kc: corner wavenumber, spectrum decays linearly for k>kc
% acf   - autocorrelation function: 
%    'gs' or 'GS' - Gaussian
%    'ex' or 'EX' - Exponential
%    'ak' or 'AK' - anisotropic vonKarman
%    'fr' or 'FR' - fractal distribution
%  Rseed - seeds for random number generators if omitted or empty
%     Rseed = sum(100*clock) is used (returned in structure spar)
%     [Rpseed Rsseed] for the phase and small random spectral part

%  OUTPUT: 2D-random field of size (Nz+1 x Nx+1) approximately standard Gaussian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4, acf = 'ak'; end
if nargin < 5, Rseed = []; end

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
        D = corr(1); kc = 0.01;
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
	   kr(i,j) = sqrt((az^2)*(kz(i)^2) + (ax^2)*(kx(j)^2)); %for PSD
       kr2(i,j) = sqrt((kz(i)^2) + (kx(j)^2)); %real wavenumber
	end
    end
end

% Power Spectral Density
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

PS = PS / max(PS(:));

if isempty(Rseed)
    rng('shuffle');
else
    rng(Rseed);
end

% Generate random complex spectrum with the prescribed PSD

[m,n] = size(kr);
xi1 = randn(m,n); xi2 = randn(m,n);
% PH = 2*pi*rand(m,n);	 
spectrum = sqrt(0.5*PS) .* (xi1 + 1i*xi2); %cartesian form
% spectrum=sqrt(PS).*exp(1i*PH);
spectrum(rmz+1,rmx+1)= sqrt(PS(rmz+1,rmx+1)) * randn; 
% Enforce Hermitian symmetry for real-valued output - Now we have 4 quadrants
U=HermitianSymmetry(spectrum,rmx,rmz);

% Inverse FFT
Y = real(ifft2(U));
Y=Y-mean(Y(:)); %zero mean
Y = Y / std(Y(:)); % to have unit variance
if mean(Y(:)) < 0; Y = (-1)*Y; end

%Resampling to get required number of points
[x_out,z_out]=meshgrid(linspace(0,N(2),nptsX),linspace(0,N(1),nptsZ));
[x,z]=meshgrid(linspace(0,N(2),2*rmx+1),linspace(0,N(1),2*rmz+1));
Y = interp2(x,z,Y,x_out,z_out,'linear');

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
