function [slip] = ConvertTruncated(slipp,slip_max)

% By Victor Hernández (victorh@hi.is). July 2025

% INPUTS:
% slipp: 2D slip distribution
% slip_max: maxmum slip from a scaling relation or an user defined value

% OUTPUT:
% slip transformed into a truncated exponential distribution (Thingbaijam and Mai, 2016) 
% using a rank transform

% Thingbaijam, K. K. S., and P. M. Mai (2016). Evidence for
% truncated exponential probability distribution of earthquake slip,
% Bull. Seismol. Soc. Am. 106, no. 4, 1802–1816, doi: 10.1785/0120150291.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    slip_mean = mean(slipp(:));
 
    % Define PDF of truncated exponential
    truncexpon_pdf = @(x, uc) (1 / uc) * exp(-x / uc) ./ (1 - exp(-slip_max / uc));
    
    % Define mean of truncated exponential using integration
    truncexpon_mean_num = @(uc) integral(@(x) x .* truncexpon_pdf(x, uc), 0, slip_max);
    ep = @(uc) slip_mean - truncexpon_mean_num(uc);
    uc = fzero(ep, [0.001, 1.5*slip_max]); 
    
    slip_flat = slipp(:);
    % Generate sorted truncated exponential samples
    samples = sort(truncexpon_sample(slip_max / uc, uc, numel(slipp)));
    % Get the sort order of the original slip values
    [~, sort_idx] = sort(slip_flat);
    % Replace the sorted values with sorted samples
    slip_flat(sort_idx) = samples;
    % Reshape back to original slip matrix size
    slip = reshape(slip_flat, size(slipp));

    function samples = truncexpon_sample(b, scale, N)
    % Truncated exponential samples from 0 to b*scale
    u = rand(N,1);
    samples = -scale.* log(1 - u .* (1 - exp(-b)));
    end

end