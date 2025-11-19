%%% Coded by Song, S. (April 2013)
% Modified by Victor Hernández (victorh@hi.is). July 2025
%%% Generate 2D distributions of source parameters
%%% using the Cholesky factorization
%%% based on 1-point and 2-point statistics

function [rup] = gen_dist_2_updated(rup)

N  = rup.nx*rup.nz;
N1 = (rup.nx1)*(rup.nz1);

str_N = sprintf('=> Number of subfaults (fine grid): %d',N);
disp(str_N), disp('  ')

str_N = sprintf('=> Number of subfaults (coarse grid): %d',N1);
disp(str_N), disp('  ')

[XX ZZ] = meshgrid(rup.lx1{1},rup.lz1{1});
[X1 Z1] = meshgrid(rup.lx{1},rup.lz{1});
    
slip_sampled   = interp2(X1,Z1,rup.slip_deterministic,XX,ZZ,'linear');
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
    rho = gen_rho_2(rup,r); %clear r

    Cm = [rho{1}{1}  rho{1}{2}  rho{1}{3}; ...
              rho{1}{2}' rho{2}{2}  rho{2}{3}; ...
              rho{1}{3}' rho{2}{3}' rho{3}{3}];  %clear rho

    %removing negative eigen values for positivity constraints
    % [eig_v, eig_d] = eig(Cm);
    % rup.eigen = diag(eig_d);
    % eigen = rup.eigen;
    % eigen(eigen<0.001) = 0.001;
    % Cm = eig_v*diag(eigen)*(eig_v)';
    % Cm = (Cm + Cm') / 2;

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
Z1 = norminv(F_map);  % Now Gaussian
% mean(Z1(:))
% std(Z1(:))

 C11 = rho{1}{1};
 C12 = rho{1}{2};
 C13 = rho{1}{3};
 C22 = rho{2}{2};
 C23 = rho{2}{3};
 C33 = rho{3}{3};

 rho12 = C12 ./ sqrt(C11 .* C22);
 rho13 = C13 ./ sqrt(C11 .* C33);

 % The C11 is from the original covariance model, however the slip was 
 %  obtained separately from a VK PSD. Recalculation of C11:
 C11_new= vk_covariance(r, rup.Hurst, rup.ax, rup.az);
 %removing negative eigen values for positivity constraints
 % [eig_v1, eig_d1] = eig(C11_new);
 % eig_d1=diag(eig_d1);
 % eig_d1(eig_d1<0.001) = 0.001;
 % C11_new = eig_v1*diag(eig_d1)*(eig_v1)';
 % C11_new = (C11_new + C11_new') / 2;
 % epsilon = 1e-8;
 % C11 = C11 + epsilon * eye(size(C11));

 %Recalculation of cross-correlations:
C12_new=rho12 .* sqrt(C11_new .* C22);
C13_new=rho13 .* sqrt(C11_new .* C33);

Cm_new=[C11_new  C12_new  C13_new;...
                 C12_new'     C22         C23;...
                 C13_new'     C23'        C33];

% clear  C11 C12 C13 C22  C23 C33 C11_new C12_new C13_new Cm rho12 rho13

 Cm_new=real(Cm_new);
 [eig_v1, eig_d1] = eig(Cm_new);
 eig_d1=diag(eig_d1);
 eig_d1(eig_d1<0.001) = 0.001;
 Cm_new = eig_v1*diag(eig_d1)*(eig_v1)';
 Cm_new = (Cm_new + Cm_new') / 2;

 C11 = Cm_new(1:N11, 1:N11);
 C12 = Cm_new(1:N11, N11+1:end);
 C21 = Cm_new(N11+1:end, 1:N11);
 C22 = Cm_new(N11+1:end, N11+1:end);

  % Compute conditional mean and covariance
 mu_cond = C21 / C11 * Z1;
 sigma_cond = C22 - C21 / C11 * C12;
 sigma_cond = (sigma_cond + sigma_cond') / 2;
 % [eig_v2, eig_d2] = eig(sigma_cond);
 % eig_d2=diag(eig_d2);


L_cond = chol(sigma_cond, 'lower'); 

randn('state',rup.seed.seed);

for iter=1:rup.num
% Sample from conditional distribution
z_cond = mu_cond + L_cond * randn(2*N1, 1);

  % s0 = randn(3*N1,1);
  % s = L*s0; %clear L

  % slip2 = s(1:N1);
  % Vr2   = s(N1+1:2*N1);
  % psv2 = s(2*N1+1:end);
  Vr1   = z_cond(1:N1);
  psv1 = z_cond(N1+1:end);

 
  Vr1   = reshape(Vr1  ,rup.nz1,rup.nx1);
  psv1  = reshape(psv1,rup.nz1,rup.nx1);
  % slip2 = reshape(slip2,rup.nz1,rup.nx1);
  % Vr2  = reshape(Vr2  ,rup.nz1,rup.nx1);
  % psv2  = reshape(psv2,rup.nz1,rup.nx1);

  % rup.slip2.dist{iter} = slip2; 
  % rup.Vr2.dist{iter}   = Vr2;    
  % rup.psv2.dist{iter} = psv2;

  rup.Vr1.dist{iter}   = Vr1;    
  rup.psv1.dist{iter} = psv1;

end  

