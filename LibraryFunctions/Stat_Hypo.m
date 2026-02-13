function [slip_pos,hypo_x,hypo_z,sampled_hypo_x,sampled_hypo_z] = Stat_Hypo(rup)
%---------------------------------------------------------
% Script to calculate probabilistic hypocenter positions based on Mai et al.(2005).
%
% INPUT
%        slip - slip values on the fault plane
%        x_loc_matrix and z_loc_matrix - matrices with coordinates of
%        subfaults (meshgrid)
%        mech - faulting mechanism ['dp' -> general dip-slip]
%                                  ['dc' -> crustal dip-slip]
%                                  ['ds' -> subduction dip-slip]
%                                  ['ss' -> strike-slip]
%                                  ['al' -> global]
%
% OUTPUT
%        slip_pos - 2D PDF for the location of the hypocenter
%        [hypo_x,hypo_z] - hypocenter coordinates with the largest probability of occurrence
%        [sampled_hypo_x,sampled_hypo_z] - sampled location from the PDF
%
% *****
%
% Author: Walter Imperatori, June 2008. 
% Modified by Victor Hernández (victorh@hi.is). July 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slip=rup.slip_deterministic;
mech=rup.mech;

% fault grid-dimensions
dim=size(slip);
dim_z=dim(1);
dim_x=dim(2);

% !samp assumed to be isotropic !
step=1;

% computes pdf for on-plane position

% --> downdip
if (mech == 'dp')
   hypz=gampdf(linspace(0,1,dim_z),7.364,0.072);
elseif (mech == 'dc')
   hypz=wblpdf(linspace(0,1,dim_z),0.692,3.394);
elseif (mech == 'ds')
   hypz=gampdf(linspace(0,1,dim_z),12.658,0.034);
elseif (mech == 'ss')
   hypz=wblpdf(linspace(0,1,dim_z),0.626,3.921);
elseif (mech == 'al')
   hypz=wblpdf(linspace(0,1,dim_z),0.612,3.353);
end

% --> along strike
hypx=normpdf(linspace(0,1,dim_x),0.5,0.23);

% calculates and scales to 1 probability based on on-plane position
hypz=hypz';
hyp_pos=hypz*hypx;
hyp_pos=hyp_pos/max(max(hyp_pos));

%------------------------------------------------------
% computes pdf based on maximum and mean slip ratio

% --> computes maximum and mean slip values
[notzero_z,notzero_x]=find(slip > 0.0);
for i=1:(max(size(notzero_z)))             %exctracts where slip is not zero
    notzero_slip(i)=slip(notzero_z(i),notzero_x(i));
end
mean_slip=mean(notzero_slip);
max_slip=max(max(slip));

% --> computes pdf for mean slip ratio
tmp_slip=slip/mean_slip;
if(mech == 'dp' | mech == 'dc' | mech == 'ds')
   hyp_mean=gampdf(tmp_slip,2.616,0.623);
elseif (mech == 'ss')
   hyp_mean=gampdf(tmp_slip,1.928,0.868);
elseif (mech == 'al')
   hyp_mean=gampdf(tmp_slip,2.210,0.748);
end

% --> computes pdf for maximum slip ratio
tmp_slip=slip/max_slip;
if(mech == 'dp' | mech == 'dc' | mech == 'ds')
   hyp_max=wblpdf(tmp_slip,0.454,1.874);
elseif (mech == 'ss')
   hyp_max=wblpdf(tmp_slip,0.446,1.551);
elseif (mech == 'al')
   hyp_max=wblpdf(tmp_slip,0.450,1.688);
end

% --> computes probability based on both previous slip ratios and scales to 1
hyp_slip=hyp_max.*hyp_mean;
hyp_slip=hyp_slip/max(max(hyp_slip));

% -------------------------------------------------------
% computes probability based on both on-plane position and slip ratios
% slip_pos=hyp_slip.*hyp_pos;  
% slip_pos=slip_pos/max(max(slip_pos));   % scales it to 1 
slip_pos=hyp_pos;

% -------------------------------------------------------
% define asperities regions (definitions from Mai et al. 2005)
[out_asp_dip,out_asp_strike]=find(slip < 0.33*max_slip);
[asp_dip,asp_strike]=find(0.33*max_slip <= slip & slip < 0.66*max_slip);
[big_asp_dip,big_asp_strike]=find(slip >= 0.66*max_slip);

% defines asperities positions on the fault plane (0=no asperity, 1=asperity, 2=big asperity)
asp_fault(1:dim_z,1:dim_x)=0;
for i=1:(max(size(asp_dip)))
    asp_fault(asp_dip(i),asp_strike(i))=1;
end
for i=1:(max(size(big_asp_dip)))
    asp_fault(big_asp_dip(i),big_asp_strike(i))=2;
end

% -------------------------------------------------------------------------------------------
% applies other constraints (paragraph "Distance between Hypocenter and Closest Asperities" of
% Mai et al. (2005)) --------> to be improved!


% --> refreshes hypocenter locations probability according to position inside/outside asperities
for i=1:(max(size(out_asp_dip)))
    slip_pos(out_asp_dip(i),out_asp_strike(i))=slip_pos(out_asp_dip(i),out_asp_strike(i))*0.48;
%    mask(out_asp_dip(i),out_asp_strike(i))=mask(out_asp_dip(i),out_asp_strike(i))*0.48;
end
for i=1:(max(size(asp_dip)))
    slip_pos(asp_dip(i),asp_strike(i))=slip_pos(asp_dip(i),asp_strike(i))*0.35;
%    mask(asp_dip(i),asp_strike(i))=mask(asp_dip(i),asp_strike(i))*0.35;
end
for i=1:(max(size(big_asp_dip)))
    slip_pos(big_asp_dip(i),big_asp_strike(i))=slip_pos(big_asp_dip(i),big_asp_strike(i))*0.16;
%    mask(big_asp_dip(i),big_asp_strike(i))=mask(big_asp_dip(i),big_asp_strike(i))*0.16;
end

slip_pos=slip_pos/max(max(slip_pos));   % scales to 1


 %% Samplig the location from the PDF

% [X, Z] = meshgrid(x, y);
X=rup.x_loc_matrix;
Z=rup.z_loc_matrix;
hypo_x = X(slip_pos == 1);
hypo_z = Z(slip_pos == 1);

slip_pos = slip_pos / sum(slip_pos(:));
pdf_flat = slip_pos(:);

% Create cumulative distribution
cdf = cumsum(pdf_flat);

% Draw a random number in [0,1)
r = rand;

% Find index where CDF exceeds r
idx = find(cdf >= r, 1, 'first');

% Convert linear index back to 2D indices
[row_idx, col_idx] = ind2sub(size(slip_pos), idx);

% Map to coordinate space
sampled_hypo_x = X(row_idx, col_idx);
sampled_hypo_z = Z(row_idx, col_idx);

end



