% By Victor Hernández (victorh@hi.is), July 2025
% Based on code by Song, S. (April 2013)
% Rupture model generator based on 1-point and 2-point statistics
% of key kinematic source parameters

function [rup] = gen_rup(rup)

[rup] = make_grid(rup);

if strcmp(rup.method,'Song')
    [rup] = gen_stats_Song(rup);
elseif strcmp(rup.method,'Savran')
    rup.p1.mu =  [mean(rup.slip_nontaper(:)) rup.user_stats.mean_Vr rup.user_stats.mean_Vmax];
    rup.p1.sig =  [std(rup.slip_nontaper(:)) rup.user_stats.sigma_Vr rup.user_stats.sigma_Vmax];
else
    disp('Error, method should Song or Savran')
    return;
end

% Step I: 2-point statistics
% generate 2D distributions of source parameters
% assuming multi-variate Gaussian with zero mean and unit std

rup = create_fields(rup);

%%%Interpolation on target grid
for iter=1:rup.num
    [X Z] = meshgrid(rup.lx{iter},rup.lz{iter});
    [X1 Z1] = meshgrid(rup.lx1{iter},rup.lz1{iter});
    rup.Vr.dist{iter}   = interp2(X1,Z1,rup.Vr1.dist{iter},X,Z,'spline');
    rup.psv.dist{iter}  = interp2(X1,Z1,rup.psv1.dist{iter},X,Z,'spline');
    rup.PT.dist{iter}  = interp2(X1,Z1,rup.PT1.dist{iter},X,Z,'spline');
end

% Step II: 1-point statistics

for iter=1:rup.num
    %%% Slip  - Slip is not modified since it is prescribed earlier
    rup.slip.dist{iter} = rup.slip_deterministic;
    rup.p1.mu1(iter,1) = mean(rup.slip.dist{iter}(:));
    rup.p1.sig1(iter,1) = std(rup.slip.dist{iter}(:));
    [rup.Mo(iter), rup.Mw(iter)] = fmomentN2(rup.slip.dist{iter},[rup.nz*rup.dz rup.nx*rup.dx]);

    %%% Rupture Velocity - Assuming a Lévy distribution
    rup.Vr.dist{iter} = rank_transform_to_Levy(rup.Vr.dist{iter},rup.user_stats.alpha,rup.user_stats.beta,rup.user_stats.gamma,rup.user_stats.delta);
    % A Gaussian distribution could be used with the next line
    % rup.Vr.dist{iter} = rup.p1.sig(2)*rup.Vr.dist{iter} + rup.p1.mu(2);
    rup.Vr.dist{iter}(rup.Vr.dist{iter}<rup.p1.min(2)) = rup.p1.min(2);
    rup.Vr.dist{iter}(rup.Vr.dist{iter}>rup.p1.max(2)) = rup.p1.max(2);
    rup.p1.mu1(iter,2) = mean(rup.Vr.dist{iter}(:));
    rup.p1.sig1(iter,2) = std(rup.Vr.dist{iter}(:));

    %%% Peak Slip Velocity - Gaussian distribution
    rup.psv.dist{iter} = rup.p1.sig(3)*rup.psv.dist{iter} + rup.p1.mu(3);
    rup.psv.dist{iter}(rup.psv.dist{iter}<rup.p1.min(3)) = rup.p1.min(3);
    rup.psv.dist{iter}(rup.psv.dist{iter}>rup.p1.max(3)) = rup.p1.max(3);
    %Limit slip/Vmax ratio to 4
    ratio=rup.slip.dist{iter}(:)./rup.psv.dist{iter}(:);
    psv_flat=rup.psv.dist{iter}(:);
    psv_flat(ratio>4)=rup.slip.dist{iter}(ratio>4)/4;
    rup.psv.dist{iter}=reshape(psv_flat,rup.nz,rup.nx);
    rup.psv.dist{iter}(rup.psv.dist{iter}<rup.p1.min(3)) = rup.p1.min(3);
    rup.psv.dist{iter}(rup.psv.dist{iter}>rup.p1.max(3)) = rup.p1.max(3);
    rup.p1.mu1(iter,3) = mean(rup.psv.dist{iter}(:));
    rup.p1.sig1(iter,3) = std(rup.psv.dist{iter}(:));

    %%% Peak Time - Computed always but not used for exponential STF
    rup.PT.dist{iter} = rup.user_stats.sigma_PT*rup.PT.dist{iter} + rup.user_stats.mean_PT;
    rup.PT.dist{iter}(rup.PT.dist{iter}<rup.p1.min(5)) = rup.p1.min(5);
    rup.PT.dist{iter}(rup.PT.dist{iter}>rup.p1.max(5)) = rup.p1.max(5);
    rup.p1.mu1(iter,5) = mean(rup.PT.dist{iter}(:));
    rup.p1.sig1(iter,5) = std(rup.PT.dist{iter}(:));

    switch rup.svf.type %Source time function STF
        case 'rec'
            rup.risT.dist{iter} = rup.slip.dist{iter}./rup.psv.dist{iter};
            rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
            rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
            rup.p1.mu1(iter,4) = mean(rup.risT.dist{iter}(:));
            rup.p1.sig1(iter,4) = std(rup.risT.dist{iter}(:));
        case 'tri'
            rup.risT.dist{iter} = rup.slip.dist{iter}./rup.psv.dist{iter}*2;
            rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
            rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
            rup.p1.mu1(iter,4) = mean(rup.risT.dist{iter}(:));
            rup.p1.sig1(iter,4) = std(rup.risT.dist{iter}(:));
        case 'pliu'
            C_pliu = 1.4*pi* rup.PT.dist{iter} + 1.2* rup.PT.dist{iter} + 0.3*pi*(1- rup.PT.dist{iter});
            rup.risT.dist{iter} = 2*rup.slip.dist{iter}*pi./rup.psv.dist{iter}./C_pliu;
            rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
            rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
            rup.p1.mu1(iter,4) = mean(rup.risT.dist{iter}(:));
            rup.p1.sig1(iter,4) = std(rup.risT.dist{iter}(:));
        case 'etinti'
            rup.risT.dist{iter} = (1.04*rup.slip.dist{iter}./(rup.PT.dist{iter}.^0.54*rup.psv.dist{iter})).^(1/0.47);
            rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
            rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
            rup.p1.mu1(iter,4) = mean(rup.risT.dist{iter}(:));
            rup.p1.sig1(iter,4) = std(rup.risT.dist{iter}(:));
        case 'exp'
            rup.risT.dist{iter} = 6.64*rup.slip.dist{iter}./(exp(1)*rup.psv.dist{iter});
            rup.risT.dist{iter}(rup.risT.dist{iter}<rup.p1.min(4)) = rup.p1.min(4);
            rup.risT.dist{iter}(rup.risT.dist{iter}>rup.p1.max(4)) = rup.p1.max(4);
            rup.p1.mu1(iter,4) = mean(rup.risT.dist{iter}(:));
            rup.p1.sig1(iter,4) = std(rup.risT.dist{iter}(:));
        otherwise
            error('The type of svf is not supported yet')
    end
end

% Step III: post-processing of generated source models
% a)Tr calculation from simulated local Vr using fast marching algorithm
disp('   ')
disp('### Step III: Post-processing   ')
disp('###   ')
vr_matrix=reshape(rup.Vs,[rup.nz,rup.nx]);
for iter=1:rup.num
    cell_hypo = [ceil(rup.dhyp/rup.sampling(2)), ceil(rup.shyp/rup.sampling(1))];
    rup.rupT.dist{iter} = fast_marching_2d(vr_matrix.*rup.Vr.dist{iter}, cell_hypo, rup.sampling(1), rup.sampling(2));
end

% b) generating slip velocity functions (SVF)
for iter=1:rup.num
    for k=1:rup.nz
        for i=1:rup.nx
            rup.svf.nt{iter}(k,i) = ceil(1.2*rup.risT.dist{iter}(k,i)/rup.svf.dt);
            [rup.svf.svf{iter}{k}{i}, rup.svf.time{iter}{k}{i}] = ...
                gen_svf(rup.risT.dist{iter}(k,i),rup.svf.dt,rup.svf.nt{iter}(k,i),rup.svf.type,rup.PT.dist{iter}(k,i));

            rup.psv2.dist{iter}(k,i) = max(rup.svf.svf{iter}{k}{i})*rup.slip.dist{iter}(k,i);
        end
    end
end

% c) generating moment rate function
Area = rup.dx*rup.dz;
for iter=1:rup.num
    tlen = ceil(max(rup.rupT.dist{iter}(:))+max(rup.risT.dist{iter}(:)))+2;
    nt = ceil(tlen/rup.svf.dt);
    mrf_time = [0:rup.svf.dt:rup.svf.dt*(nt-1)];
    mrf = zeros(rup.nx*rup.nz,nt);

    mu=reshape(rup.mu_vec,rup.nz,rup.nx);
    ind = 0;
    for k=1:rup.nz
        for i=1:rup.nx
            ind = ind + 1;
            i1 = round(rup.rupT.dist{iter}(k,i)/rup.svf.dt)+2;
            i2 = i1 + length(rup.svf.svf{iter}{k}{i}) - 1;
            mrf(ind,i1:i2) = mu(k,i)*Area*rup.svf.svf{iter}{k}{i}*rup.slip.dist{iter}(k,i);
        end
    end
    rup.mrf.mrf{iter}  = sum(mrf);
    rup.mrf.time{iter} = mrf_time;
    [~,rup.mrf.fas{iter},rup.mrf.freq{iter}]= Compute_Fourier(rup.mrf.mrf{iter},rup.svf.dt);
end

disp('#### Finished ####');




