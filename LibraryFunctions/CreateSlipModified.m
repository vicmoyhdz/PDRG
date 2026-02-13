function [slip, slip_taper,slipp,radial_spectra,Deterministic_resampled] = CreateSlipModified(Deterministic, wavelength_cut, slip_max,dim, Mw, mu, samp, cor_len, acf_type, taper_width, taper_function,seed)

% By Victor Hernández (victorh@hi.is). July 2025
% This function merges a given slip distribution (Deterministic) with an stochastic 
% distribution created using the spectral synthesis method (SpecSyn2D_Hybrid).
% For merging uses the filter by Graves and Pitarka (2010). Then, it
% modifies slip values to follow a truncated exponential
% distribution, and later scales the field to the required seismic moment.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 11
        taper_width = [1000, 1000, 1000];
    end
    if nargin < 12
        taper_function = 'hn';
    end

    [F,radial_spectra,Deterministic_resampled] = SpecSyn2D_Hybrid(Deterministic, wavelength_cut,dim, samp, cor_len, acf_type,seed);

    W = dim(1);  % Width in km
    L = dim(2);  % Length in km

    % Lift field above zero
    Y = F - min(F(:));

    % Rescale slip to match desired seismic moment
    [mo, ~] = fmomentN(Y, [W, L], mu);
    Mo = 10^(1.5*Mw + 9.05);
    slipp = Y * (Mo / mo);


    % -- Apply taper --
    ntx = round(taper_width(1) / samp(2));
    ntt = round(taper_width(2) / samp(1));
    ntb = round(taper_width(3) / samp(1));
    slipp = TaperSlip(slipp, [ntx, ntt, ntb], taper_function);

    % slip_mean = mean(slipp(:));
    % slip_max = 10^(0.95 * log10(slip_mean) + 0.52);  % umax from Figure 9 Thingbaijam and Mai

   [slip] = ConvertTruncated(slipp,slip_max);
   % slip=slipp;

    % -- Apply taper --
    ntx = round(taper_width(1) / samp(2));
    ntt = round(taper_width(2) / samp(1));
    ntb = round(taper_width(3) / samp(1));
    slip_taper = TaperSlip(slip, [ntx, ntt, ntb], taper_function);

    % -- Rescale after taper to restore moment --
    [mo, ~] = fmomentN(slip_taper, [W, L],mu);
    Mo = 10^(1.5*Mw + 9.05);
    slip_taper = slip_taper * (Mo / mo);
end

function [Mo,Mw] = fmomentN(ina,dim,mu)

    W = dim(1); L = dim(2);
    [Nz,Nx]=size(ina);
    dx=L/Nx; dz=W/Nz;
    A = dx*dz;			% subfault surface [m]
    
    % mu = 3.3*1e10;			% shear modulus = rigidity [N/m^2]
    % s = mean(ina,'all');		% average slip over fault surface [m]
    % A = L*W;			% fault surface [m]  
    % Mo = mu*s*A;			% moment [Nm]

    Mo=sum(A*mu.*ina(:));
    Mw = (2/3)*(log10(Mo) - 9.05);	% moment magnitude
end


function S = TaperSlip(S, N, window)
% TaperSlip tapers the slip amplitudes at fault boundaries
% 
% INPUT:
%   S       - Original 2D slip grid
%   N       - [left/right, top, bottom] taper widths (in points)
%   window  - 'hn' (Hanning), 'kw' (Kaiser), 'tr' (Triangular), or a custom vector
%
% OUTPUT:
%   S       - Tapered slip grid
%
% Author: Martin Mai
% Original: 1998-1999, Stanford

if nargin < 2
    N = [1 0 1];
end
if nargin < 3
    window = 'hn';
end

[m, n] = size(S);
bell = ones(m, n);  % bell-shaped taper window

% --------------------------
% Generate taper windows
% --------------------------
if ischar(window)
    switch window
        case 'hn'
            taperS = hanning(2 * N(1) + 3); taperS = taperS(2:end-1);
            taperT = hanning(2 * N(2) + 3); taperT = taperT(2:end-1);
            taperB = hanning(2 * N(3) + 3); taperB = taperB(2:end-1);
        case 'kw'
            beta = 6;
            taperS = kaiser(2 * N(1) + 1, beta);
            taperT = kaiser(2 * N(2) + 1, beta);
            taperB = kaiser(2 * N(3) + 1, beta);
        case 'tr'
            taperS = triang(2 * N(1) + 1);
            taperT = triang(2 * N(2) + 1);
            taperB = triang(2 * N(3) + 1);
        otherwise
            error('Unknown window type');
    end

    % Extract one-sided tapers
    winS = taperS(N(1)+2:end);
    winT = taperT(1:N(2));
    winB = taperB(end-N(3)+1:end);

else
    % custom taper vector case (must be scalar multipliers)
    for i = 1:length(N)
        if N(i) == 0
            window(i) = 1;
        end
    end
    winS = window(1);
    winT = window(2);
    winB = window(3);
end

% --------------------------
% Apply side taper (left/right)
% --------------------------
if length(winS) > 1
    ls = length(winS);
    for row = 1:m
        bell(row, 1:ls) = bell(row, 1:ls) .* flip(winS(:))';
        bell(row, end-ls+1:end) = bell(row, end-ls+1:end) .* winS(:)';
    end
elseif isscalar(winS)
    bell = bell * winS;
end

% --------------------------
% Apply top/bottom taper
% --------------------------
if length(winT) > 1
    lt = length(winT);
    for col = 1:n
        bell(1:lt, col) = bell(1:lt, col) .* winT(:);
    end
end

if length(winB) > 1
    lb = length(winB);
    for col = 1:n
        bell(end-lb+1:end, col) = bell(end-lb+1:end, col) .* winB(:);
    end
end

% Final tapered slip
S = S .* bell;
end
