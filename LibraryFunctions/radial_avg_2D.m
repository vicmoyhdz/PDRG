function [krad, Pkrad] = radial_avg_2D(KR, S2D)
    nBins = length(KR)/2;
    % KR = sqrt(kx.^2 + ky.^2);        % circular radial wavenumber
    kmax = max(KR(:));
    edges = linspace(0, kmax, nBins+1);
    krad = 0.5*(edges(1:end-1)+edges(2:end));
    Pkrad = nan(size(krad));
    for i = 1:numel(krad)
        mask = (KR >= edges(i)) & (KR < edges(i+1));
        if any(mask(:))
            Pkrad(i) = mean(S2D(mask));
        end
    end
end