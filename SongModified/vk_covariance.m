function C = vk_covariance(r, H, ax, az, sigma2)
    % r contains x and z distance lags
    % H: Hurst exponent
    % ax, az: correlation lengths
    % sigma2: sill (variance)

    if nargin < 5
        sigma2 = 1;
    end

    d = sqrt((r.x / ax).^2 + (r.z / az).^2);

    vk_cov = @(h, H) sigma2 * (2^(1 - H) / gamma(H)) * ...
    (h+1e-6).^H .* besselk(H,h+1e-6);

   C=vk_cov(d,H);
   C(d < 1e-6) = sigma2;
