%%%By Victor Hernández (victorh@hi.is), July 2025
%%% Modified from code by Song, S. (April 2013),
%%% Generate 2D distributions of source parameters
%%% using the Cholesky factorization conditioned on slip
%%% based on 1-point and 2-point statistics

function [rup] = create_fields(rup)

N  = rup.nx*rup.nz;
N1 = (rup.nx1)*(rup.nz1);

[XX ZZ] = meshgrid(rup.lx1{1},rup.lz1{1});
[X1 Z1] = meshgrid(rup.lx{1},rup.lz{1});

slip_sampled   = interp2(X1,Z1,rup.slip_nontaper,XX,ZZ,'linear');
slip_sampled(1,:)=slip_sampled(2,:);
slip_sampled(end,:)=slip_sampled(end-1,:);
slip_sampled(:,1)=slip_sampled(:,2);
slip_sampled(:,end)=slip_sampled(:,end-1);

X = single(XX(:));
Z = single(ZZ(:));

Z1 = repmat(Z ,1,N1); clear Z
r.z = (Z1' - Z1);    clear Z1

X1 = repmat(X ,1,N1); clear X
r.x = (X1' - X1);    %clear X1

N11= length(r.x);

if strcmp(rup.method,'Song')
    rho = gen_rho(rup,r); %clear r

    Cm = [rho{1}{1}  rho{1}{2}  rho{1}{3}; ...
        rho{1}{2}' rho{2}{2}  rho{2}{3}; ...
        rho{1}{3}' rho{2}{3}' rho{3}{3}];  %clear rho

    %removing negative eigen values for positivity constraints
    [eig_v, eig_d] = eig(Cm);
    rup.eigen = diag(eig_d);
    eigen = rup.eigen;
    eigen(eigen<0.001) = 0.001;
    Cm = eig_v*diag(eigen)*(eig_v)';
    Cm = (Cm + Cm') / 2;

else
    [rho] = covariance_coregionalization(rup,r); %Savran model
    Cm = [rho{1}{1}  rho{1}{2}  rho{1}{3}; ...
        rho{1}{2}' rho{2}{2}  rho{2}{3}; ...
        rho{1}{3}' rho{2}{3}' rho{3}{3}];  %clear rho
end

% L = chol(Cm,'lower');   clear Cm

%Convert given slip to gaussian
f_vec = slip_sampled(:);
[f_sorted, sort_idx] = sort(f_vec);
nn = length(f_sorted);
F_empirical = ((1:nn)' - 0.5) / nn;  % Empirical CDF values

% Step 2: Assign CDF values to original data
F_map = zeros(nn,1);
F_map(sort_idx) = F_empirical;

% Step 3: Transform to Gaussian using inverse normal CDF
epsilon = 1e-10;
F_map = min(max(F_map, epsilon), 1 - epsilon);  % Clamp
Z1_gaussian = norminv(F_map);  % Now Gaussian

C11 = Cm(1:N11, 1:N11);
C12 = Cm(1:N11, N11+1:end);
C21 = Cm(N11+1:end, 1:N11);
C22 = Cm(N11+1:end, N11+1:end);

% Step 4: Compute conditional mean and covariance
mu_cond = C21 / C11 * Z1_gaussian;
sigma_cond = C22 - C21 / C11 * C12;
sigma_cond = (sigma_cond + sigma_cond') / 2;
% [eig_v2, eig_d2] = eig(sigma_cond); %checking positive definitiveness
% eig_d2=diag(eig_d2);


L_cond = chol(sigma_cond, 'lower');

%Correlation matrix for peak time
C_PT= rho{3}{3};
L_PT= chol(C_PT, 'lower');

rng(rup.seed.seed);

Mean_Vmax=rup.user_stats.mean_Vmax*0.7;
rup.iterations=0;

while (Mean_Vmax<rup.user_stats.mean_Vmax*0.95) || (Mean_Vmax>rup.user_stats.mean_Vmax*1.05)

    for iter=1:rup.num
        % Sample from conditional distribution
        z_cond = mu_cond + L_cond * randn(2*N1, 1);

        Vr1   = z_cond(1:N1);
        psv1 = z_cond(N1+1:end);

        % Create distribution of peak time
        z_PT = randn(N1, 1);
        PT_std = L_PT * z_PT;
        PT1=rup.user_stats.rho_PT * psv1 + sqrt(1 - rup.user_stats.rho_PT^2) * PT_std;

        Vr1   = reshape(Vr1  ,rup.nz1,rup.nx1);
        psv1  = reshape(psv1,rup.nz1,rup.nx1);
        PT1  = reshape(PT1,rup.nz1,rup.nx1);

        rup.Vr1.dist{iter}   = Vr1;
        rup.psv1.dist{iter} = psv1;
        rup.PT1.dist{iter} = PT1;
    end

    Vmax=rup.p1.sig(3)*rup.psv1.dist{1} + rup.p1.mu(3);
    Mean_Vmax=mean(Vmax(:));
    rup.iterations=rup.iterations+1;
end
end

