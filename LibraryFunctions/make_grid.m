
function [rup] = make_grid(rup)

% Basic information for an event

rup.target_Mo = 10.^(1.5*rup.target_Mw + 9.05);

rup.dx   = rup.sampling(2);      % grid size (along-strike) in km
rup.dz   = rup.sampling(1);      % grid size (along-dip) in km

rup.dz1 =rup.dx1;    % grid size for sub-calculation
rup.svf.dt   = 0.005;
%%% End of inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rup.nx = ceil(rup.L/rup.dx);
rup.nz = ceil(rup.W/rup.dz);

rup.nx1 = ceil(rup.L/rup.dx1);
rup.nz1 = ceil(rup.W/rup.dz1);

for inum=1:rup.num
  rup.lx{inum} = [0:rup.dx:rup.dx*(rup.nx-1)];
  rup.lz{inum} = [0:rup.dz:rup.dz*(rup.nz-1)];
  
  rup.lx1{inum} = [0:rup.dx1:rup.dx1*(rup.nx1-1)];
  rup.lz1{inum} = [0:rup.dz1:rup.dz1*(rup.nz1-1)];

  rup.lx1{inum} = [(rup.lx1{inum}(1)-rup.dx1) rup.lx1{inum} (rup.lx1{inum}(end)+rup.dx1)];
  rup.lz1{inum} = [(rup.lz1{inum}(1)-rup.dz1) rup.lz1{inum} (rup.lz1{inum}(end)+rup.dz1)];
  
  [XX, ZZ] = meshgrid(rup.lx{inum},rup.lz{inum});
  rup.dis{inum} = sqrt(XX.^2 + ZZ.^2);  % distance from the nucleation point on the fault plane

end

  rup.nx1 = rup.nx1 + 2;
  rup.nz1 = rup.nz1 + 2;
  


