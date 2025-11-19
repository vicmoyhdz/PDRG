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