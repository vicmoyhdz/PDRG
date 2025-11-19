function [Roughness_subfault,normal_vectors] = fault_roughness_frac(ISD,alpha,strike_length,dip_length,dx_subfault, dy_subfault)

% By Victor Hernández (vmh5@hi.is). July 2025
% Generates fault roughness: phase is extracted from ISD field whilst the spectral amplitude
% is extracted from a fractal PSD using the spectral sythesis method (Ogorodnikov). This
% approach was proposed by Aquib et al. (2025).

% Shi, Z., and S. M. Day (2013), Rupture dynamics and ground motion from 3-D rough-fault 
% simulations, J. Geophys. Res. Solid Earth, 118, 1122–1141, doi:10.1002/jgrb.50094.

% Ogorodnikov, V.A. and Prigarin, S.M. (1996), Numerical Modelling of Random Processes 
% and Fields: Algorithms and Applications. De Gruyter.

% INPUTS:
% ISD: integrated stress field. Must cover the same are of the fault
% alpha: roughness parameter. Power and Tullis (1991) provide estimates of α in the 
% 10−3 to 10−2 range for natural faults.
% Dimensions of the fault and subfaults

% OUTPUT:
% Roughness (in m) at the center of each subfault (out-of-plane perturbation). 
% normal_vectors at the center of each subfault in the fault's coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up

H=1; %Hurst exponent ->To have a k^-2 decay
D=2; %Fractal dimension

Nx=2*strike_length/dx_subfault+1; Ny=2*dip_length/dy_subfault+1; 
% The number of points is duplicated so that the gradient can be computed
% at the center of each subfault (to compute normals)
dx=strike_length/(Nx-1); dy=dip_length/(Ny-1);
[x,y] = meshgrid((0:Nx-1) * dx, (0:Ny-1) * dy);
[x1,y1] = meshgrid((0:Nx/2-1) * dx*2+dx, (0:Ny/2-1) * dy*2+dy);

% Frequency grid
dkx = 1/ (Nx * dx);
dky = 1 / (Ny * dy);
kx = ifftshift((-floor(Nx/2):ceil(Nx/2)-1) * 2*pi*dkx);
ky = ifftshift((-floor(Ny/2):ceil(Ny/2)-1) * 2*pi*dky);
[KX, KY] = meshgrid(kx, ky);
kr =  sqrt(KY.^2 + KX.^2); %radial wavenumber
% Avoid zero to prevent division by zero in 1/k
kr(kr == 0) = 1e-12;

% --- Define PSD: 
PSD = (kr).^(-2*(H+D/2)) ;

% Convert wavelength cutoffs to wavenumber cutoffs
lambda_min = dx_subfault;    % shortest wavelength
lambda_max = strike_length;     % longest wavelength
k_max = 2*pi / lambda_min;
k_min = 2*pi  / lambda_max;

%% Generate random complex spectrum with filtered PSD and phase from given ISD field
 
%Interpolate ISD filed
ISD_sampled=interp2(x1,y1,ISD,x,y);
ISD_sampled(1,:)=ISD_sampled(2,:);
ISD_sampled(end,:)=ISD_sampled(end-1,:);
ISD_sampled(:,1)=ISD_sampled(:,2);
ISD_sampled(:,end)=ISD_sampled(:,end-1);
ISD_fft = fft2(ISD_sampled);
phase_ISD = angle(ISD_fft);

xi1 = randn(Ny,Nx); xi2 = randn(Ny,Nx);
amp = sqrt(0.5*PSD) .* sqrt(xi1.^2 + xi2.^2);   % Rayleigh-distributed amplitude
% --- Apply band-pass filter ---
band_mask = (kr >= k_min) & (kr <= k_max);
amp_filtered=amp.*band_mask;

spectrum=amp_filtered .* exp(1i * phase_ISD); %polar form

% Enforce Hermitian symmetry for real-valued output
for i = 1:Ny
    for j = 1:Nx
        i_sym = mod(Ny - i + 1, Ny) + 1;
        j_sym = mod(Nx- j + 1, Nx) + 1;
        spectrum(i, j) = conj(spectrum(i_sym, j_sym));
    end
end

% Inverse FFT to spatial domain
X = real(ifft2(spectrum));

%% Scalling of roughness
Roughness=X(1:Ny,1:Nx);
Roughness = Roughness - mean(Roughness(:));

%Scale to have requested alpha (std)
alpha_calc = std(Roughness(:))/strike_length;
Roughness=alpha*Roughness./alpha_calc;

%% Normal vectors

%Compute gradientes at the centers of the subfaults
dZddip=(Roughness(3:2:end,2:2:end)-Roughness(1:2:end-2,2:2:end))/(2*dy);
dZdstrike=(Roughness(2:2:end,3:2:end)-Roughness(2:2:end,1:2:end-2))/(2*dx);

%Compute normal vectors at the centers of the subfaults
Normx = -dZdstrike;
Normy = -dZddip;
Normz = ones(size(dZdstrike));

% Normalize
normLength = sqrt(Normx.^2 + Normy.^2 + Normz.^2);
Normx = Normx ./ normLength;
Normy = Normy ./ normLength;
Normz = Normz ./ normLength;
normal_vectors=[Normx(:),Normy(:),Normz(:)];

% Resample roughness at subfaults
 Roughness_subfault=Roughness(2:2:end,2:2:end);

end