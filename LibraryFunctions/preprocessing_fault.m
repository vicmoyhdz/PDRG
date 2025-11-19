function [rup] = preprocessing_fault(rup)

% By Victor Hernández (victorh@hi.is). July 2025
% Function to discretize the fault and compute global coordinates of the subfaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reading Mechanical Properties file
 fid = fopen(rup.mech_file,'r');
 Num_Mat = str2num(fgetl(fid));
 mec_prop = zeros(Num_Mat(1),6);
 for i = 1 : Num_Mat(1)
     mec_prop(i,:) = str2num(fgetl(fid));
 end
 fclose(fid);

dx=rup.sampling(2); dz=rup.sampling(1);
Nx=rup.L/rup.sampling(2); Nz=rup.W/rup.sampling(1);
[x,z] = meshgrid((0:Nx-1) * dx, (0:Nz-1) * dz); 
rup.x_loc_matrix=x; rup.z_loc_matrix=z;

 mec_prop(:,1:4) = mec_prop(:,1:4)*1000;
 lay_bot_surf_z = rup.z_ref-cumsum(mec_prop(:,4)); 
 x_ref=rup.sampling(2)/2+rup.x_loc_matrix;
 y_ref=rup.sampling(1)/2+rup.z_loc_matrix;
 x_loc = x_ref;
 y_loc = -y_ref*cosd(rup.dip);
 z_loc = -y_ref*sind(rup.dip);
 % Rotating Local coordinate to Global coordinate system
 [x_glo,y_glo] = find_rotated_coords(x_loc,y_loc,-(rup.stk -90));
 % Translating the Origin to Global reference Coordnate System
 x_glo = x_glo + rup.ref(1);
 y_glo = y_glo + rup.ref(2);
 z_glo = z_loc + rup.ref(3);
 rup.glo_coor = [x_glo(:) y_glo(:) z_glo(:)];
 % %Hypo coordinates
 % HypLoc = [rup.shyp -rup.dhyp*cosd(rup.dip) -rup.dhyp*sind(rup.dip)];    
 % [HypGlo(1),HypGlo(2)] = find_rotated_coords(HypLoc(1),HypLoc(2),-(rup.stk-90));
 % HypGlo(1) = HypGlo(1) + rup.ref(1);
 % HypGlo(2) = HypGlo(2) + rup.ref(2);
 % HypGlo(3) = HypLoc(3) + rup.ref(3);
 % rup.hypo_coor = HypGlo;


 for ipt = 1:length(x(:))
         indx1 = find(lay_bot_surf_z <= z_glo(ipt));
         if isempty(indx1)
             lay_num(ipt) = length(lay_bo_surf_z);
         else
            lay_num(ipt) = indx1(1);
         end
  end
  mu = mec_prop(:,3).*mec_prop(:,2).*mec_prop(:,2);
  
  rup.Rho = mec_prop(lay_num,3);
  rup.Vs = mec_prop(lay_num,2);
  rup.Vp = mec_prop(lay_num,1);
  rup.mu_vec = mu(lay_num,1);

end