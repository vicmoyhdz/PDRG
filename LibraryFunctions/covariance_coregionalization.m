function [rho] = covariance_coregionalization(rup,r)

% By Victor Hernández (victorh@hi.is). July 2025
% Function implementing the Covariance LMC model
% by Savran and Olsen (2020)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  gaussian_cov = @(h, a) exp(-(3*h) / a);

    B1=[0.0282 0.0164 0.002;
           0.0164 0.6917 0.0631;
           0.0020 0.0631 0.0403];

     B2=[0.9718 0.1504 0.8100;
            0.1504 0.3083 0.1841;
            0.8100 0.1841 0.9597];

    D = sqrt(r.x.^2+r.z.^2);  % N1xN1 distance matrix
    C1 = gaussian_cov(D, 250);
    C2 = gaussian_cov(D, 5000);  % N1xN1 each

    P = length(rup.p1.mu);  % number of fields

    for k = 1:P
        for i = 1:P
              rho{k}{i}  = B1(k,i)*C1 + B2(k,i)*C2;
        end
    end
    
end