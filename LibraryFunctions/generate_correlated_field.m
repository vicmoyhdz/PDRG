function [xi2, corr_coeff] = generate_correlated_field(S, eta,std_target,mean_target)

% By Victor Hernández (vmh5@hi.is). July 2025
% Generate_correlated_field - Generates a random field xi2 with:
%   - approximately the same PSD as the given field S
%   - correlation coefficient ≈ eta with S
%
% Inputs:
%   S   - given real-valued 2D field
%   eta - desired correlation coefficient (between -1 and 1)
%   std_target - Std of the output field
%   mean_target - Mean of the output field
%
% Outputs:
%   xi2        - generated random field with correlation ≈ eta and same PSD as S
%   corr_coeff - actual computed correlation coefficient
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% Check input
assert(eta >= -1 && eta <= 1, 'eta must be between -1 and 1');
assert(isreal(S), 'Input field S must be real-valued');

% Get size
[Nx, Ny] = size(S);

% Compute PSD of S
S_hat = fft2(S);
PSD_S = abs(S_hat).^2;

% Generate a random complex white noise field with uniform phase
E1 = randn(Nx, Ny) + 1i * randn(Nx, Ny);

% Impose same PSD
E1_scaled = sqrt(PSD_S/2) .* E1;
eps1 = real(ifft2(E1_scaled));

% Ensure zero-mean, unit variance
S_std = (S - mean(S(:))) / std(S(:));
eps1_std = (eps1 - mean(eps1(:))) / std(eps1(:));

eps1_proj = eps1_std - (dot(S_std(:), eps1_std(:)) / norm(S_std(:))^2) * S_std;

% Normalize again
eps1_proj = (eps1_proj - mean(eps1_proj(:))) / std(eps1_proj(:));

% Combine
xi2_std = eta * S_std + sqrt(1 - eta^2) * eps1_proj;

% Restore mean and std of original S
xi2 = xi2_std * std_target + mean_target;

% Compute actual correlation coefficient (for checking)
corr_coeff = corr(S(:), xi2(:));

end
