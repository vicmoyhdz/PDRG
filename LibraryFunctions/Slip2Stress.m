function [sigmaS,sigmaD] = Slip2Stress(slip,rake,samp,rig,lam,sfac)
%
%   [sigmaS,sigmaD,EsS,EsD] = Slip2Stress9(slipS,slipD,samp,rig,lam,sfac);
%
%   This function calculates the static stress drop for a given
%   slip distribution using Andrews (1980) concept of a static
%   stiffness function that involves a 2D-Fourier Transform of
%   the slip on the fault; it also returns an estimate for the
%   the static self-energy on the fault, Es
%
%   Input : slipS - strike slip distribution in (m)
%           slipD - dip slip distribution in (m)
%           samp  - sampling in x,z direction [dz dx] in (m)
%                   default: [1 1]
%           rig   - rigidity in (N/m^2),
%                   default: 3.3e10
%           lam   - lambda (Lame's constant) in (N/m^2)
%                   default: lam = rig
%           sfac  - scaling-factor for stress; nominally sfac = 0.5,
%                   but that sometimes leads to very large stress drops
%                   default: 0.5
%
%   Output: sigmaS - along strike componnent of static stress drop in [MPa]
%           sigmaD - down dip component of static stress drop in [MPa]
%           EsS, EsD - static self-energy (after Andrews, 1980) [in J/m^2]
%                      calculated from strike slip and dip slip separately
%
%   Original version Slip2Stress.m by
%   Martin Mai, 09/05/2000 (mmai@pangea.stanford.edu)
%   -------------------------------------------------
%   Slip2Stress by Johannes Ripperger, 2004

%%% Argument check and default values
if nargin<2; error('Insufficient Input Arguments.'); end
if nargin<3 | isempty(samp);  samp = [1 1];   end
if nargin<4 | isempty(rig);   rig = 3.3*1e10; end
if nargin<5 | isempty(lam);   lam = rig;      end
if nargin<6 | isempty(sfac);  sfac = 0.5;     end

% %%% Verify same size of dip slip and strike slip
% if (size(slipS,1)~=size(slipD,1)) | (size(slipS,2)~=size(slipD,2))
%   error('Array size mismatch between slip distributions!');
% end

% Replicate the first and last rows and columns (10 times)
slip_s1 = [repmat(slip(1, :), 10, 1); slip; repmat(slip(end, :), 10, 1)];
slip_s2 = [repmat(slip_s1(:, 1), 1, 10), slip_s1, repmat(slip_s1(:, end), 1, 10)];
% slip_s2=slip;

slipS=slip_s2*cosd(rake);
slipD=slip_s2*sind(rake);

%%% Calculate fault dimensions from slip array size and sampling
dim = (size(slip_s2)-1).*samp; % in m

%%% Determine FFT grid size
N = 2^(nextpow2(max(size(slip_s2))));  

%%% Get the 2D Fourier amplitude spectra
FT2DS = samp(1)*samp(2)*fft2(slipS,N,N);
FT2DD = samp(1)*samp(2)*fft2(slipD,N,N);
PHS   = angle(FT2DS);
PHD   = angle(FT2DD);
AMS   = fftshift(abs(FT2DS));
AMD   = fftshift(abs(FT2DD));

%%% Determine wavenumber vectors fx and fz
knyz = 1/(2*samp(1));  % Nyquist wavenumbers
knyx = 1/(2*samp(2));  %
[nfz,nfx] = size(AMS);  % Retrieve the actual FFT grid size
if mod(nfx,2)==0; % nfx even
  fx=linspace(-knyx,knyx,nfx+1);
  fx(nfx+1)=[]; % Value for knyx is NOT included for even nfx
else % nfx uneven
  fx=linspace(-knyx,knyx,nfx);
end
if mod(nfz,2)==0; % nfz even
  fz=linspace(-knyz,knyz,nfz+1);
  fz(nfz+1)=[]; % Value for knyz is NOT included for even nfz
else % nfz uneven
  fz=linspace(-knyz,knyz,nfz);
end

%%% Set-up grid to calculate stress drop
[kx,kz] = meshgrid(fx,fz);
k = sqrt(kx.^2 + kz.^2);

%%% Find indices of zero wavenumber and exchange k there
%%% with a dummy value to avoid division by zero. (K is
%%% set to zero at this point later, anyway).
[zz,zx]=find(k == 0);
k(zz,zx)=eps;

%%% Static stiffness function
% As originally used by Martin, assuming lam=rig:
%   KS = (sfac*rig./k).*((4/3).*kx.^2 + kz.^2);
%   KD = (sfac*rig./k).*(kx.^2 + (4/3).*kz.^2);
% As in Andrews, 1980 JGR, eq. (20):
KS = (sfac*rig./k).*((2*(lam+rig)/(lam+2*rig)).*kx.^2 + kz.^2);
KD = (sfac*rig./k).*(kx.^2 + (2*(lam+rig)/(lam+2*rig)).*kz.^2);
KSD = (sfac*rig./k).*((2*(lam+rig)/(lam+2*rig))-1).*kx.*kz;

%%% Set K to zero for kx=kz=k=0
KS(zz,zx) = 0;
KD(zz,zx) = 0;
KSD(zz,zx) = 0;

%%% Static stress drop on the fault in FT-domain
%%% (Andrews 1980, eq. (19))
SDFS = AMS.*KS;
SDFD = AMD.*KD;
SDFSD = AMS.*KSD;
SDFDS = AMD.*KSD;

%%% Static stress drop on fault in space domain (in Pa)
i = sqrt(-1);
FS = ifftshift(SDFS).*exp(i*PHS);
FD = ifftshift(SDFD).*exp(i*PHD);
FSD = ifftshift(SDFSD).*exp(i*PHS);
FDS = ifftshift(SDFDS).*exp(i*PHD);
sigmas = real(ifft2(FS))+real(ifft2(FDS));
sigmad = real(ifft2(FD))+real(ifft2(FSD));

%%% Cut to original size
sigmaS = sigmas(1:size(slipS,1),1:size(slipS,2));
sigmaD = sigmad(1:size(slipS,1),1:size(slipS,2));

%%% Correction factor
%%% (for the MATLAB-way of defining the FFT?)
sigmaS = 2*pi*sigmaS;
sigmaD = 2*pi*sigmaD;

%%% Account for spatial sampling
sigmaS = sigmaS/(samp(1)*samp(2));
sigmaD = sigmaD/(samp(1)*samp(2));

%%% From Pa to MPa
sigmaS = sigmaS./1e6;
sigmaD = sigmaD./1e6;
sigmaS=sigmaS(11:end-10,11:end-10);
sigmaD=sigmaD(11:end-10,11:end-10);
% sigmaS(1,:)=sigmaS(2,:);
% sigmaS(end,:)=sigmaS(end-1,:);
% sigmaS(:,1)=sigmaS(:,2);
% sigmaS(:,end)=sigmaS(:,end-1);
% sigmaD(1,:)=sigmaD(2,:);
% sigmaD(end,:)=sigmaD(end-1,:);
% sigmaD(:,1)=sigmaD(:,2);
% sigmaD(:,end)=sigmaD(:,end-1);

% %%% Compute static self-energy, using method by Andrews (1980):
% %%% 'Fault impedance and earthquake energy in the fourier
% %%%  transform domain', BSSA, Vol 70, No.5, pp.1683-1698
% %%% Wavenumber sampling intervals:
% dk = [fx(2)-fx(1); fz(2)-fz(1)]; % [1/m]
% %%% Equation (52):
% E52S = 0.5*sum(sum(KS.*(AMS.^2)))*dk(1)*dk(2); % in [J]
% E52D = 0.5*sum(sum(KD.*(AMD.^2)))*dk(1)*dk(2); % in [J]
% %%% Correction factor (FFT ?)
% E52S=E52S*2*pi;
% E52D=E52D*2*pi;
% %%% We want it in J/m^2 for comparison between different events:
% EsS = E52S/(dim(1)*dim(2)*1e6); % in [J/m^2]
% EsD = E52D/(dim(1)*dim(2)*1e6); % in [J/m^2]



