function [slip] = ConvertTruncatedCauchy(slipp,slip_max)

% By Victor Hernández (victorh@hi.is). July 2025

% INPUTS:
% slipp: 2D slip distribution
% slip_max: maxmum slip from a scaling relation or an user defined value

% OUTPUT:
% slip transformed into a truncated Cauchy distribution (Lavallée et al., 2006)
% using a rank transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    slip_mean = mean(slipp(:));
    D0=0.5*slip_mean;
 
    % Wrapper function for fzero (k is the only free parameter)
    objective = @(k) cauchy_mean(k, D0, slip_max) - slip_mean;

    k_opt = fzero(objective, [0.01, 100]);

    computed_mean = cauchy_mean(k_opt,D0,slip_max);
    
    slip_flat = slipp(:);
    slip_flat = rank_transform_to_truncated_cauchy(slip_flat, D0, k_opt, slip_max);

    slip = reshape(slip_flat, size(slipp));

end

    function mu = cauchy_mean(k, D0, Dmax)
    % Unnormalized Cauchy PDF
    f = @(D) 1 ./ (1 + ((D - D0) ./ k).^2);

    % Normalization constant
    C = 1 / integral(f, 0, Dmax);

    % Normalized PDF
    pdf = @(D) C * f(D);

    % Compute the mean numerically
    mu = integral(@(D) D .* pdf(D), 0, Dmax);
end

   function y = rank_transform_to_truncated_cauchy(x, D0, k, Dmax)
    % Rank-transform vector x to match a truncated Cauchy distribution
    % using f(D) = C / (1 + ((D - D0)/k)^2), D in [0, Dmax]

    N = numel(x);

    % Step 1: Sort x and get rank order
    [~, sort_idx] = sort(x);           % Sort x
    [~, ranks] = sort(sort_idx);       % Ranks of each element in x

    % Step 2: Generate target quantiles evenly spaced in [0, 1]
    % Avoid 0 and 1 to prevent infinities in tails
    p = ((1:N)' - 0.5) / N;  % Empirical quantiles

    % Step 3: Build CDF from your truncated Cauchy
    D_vals = linspace(0, Dmax, 10000);  % Support
    unnorm_pdf = 1 ./ (1 + ((D_vals - D0)./k).^2);
    norm_const = trapz(D_vals, unnorm_pdf);
    pdf_vals = unnorm_pdf / norm_const;
    cdf_vals = cumtrapz(D_vals, pdf_vals);
    cdf_vals = cdf_vals / max(cdf_vals);  % Normalize to [0, 1]

    % Step 4: Build inverse CDF numerically
    inv_cdf = @(u) interp1(cdf_vals, D_vals, u, 'linear', 'extrap');

    % Step 5: Map quantiles to Cauchy values
    cauchy_vals = inv_cdf(p);

    % Step 6: Assign values based on rank
    y = cauchy_vals(ranks);
end