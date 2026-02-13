function [Mo,Mw] = fmomentN2(in,dim,mu)
%
% function [Mo,Mw] = fmomentN(in,dim) calculates the 
% moment of an event for a given slip distribution (in cm)
% if only IN is given, DIM = size(IN)
%
% INPUT:  in 	- 2D-array that constains slip values
%	      dim	- source dimensions [W L]
%
% OUTPUT: Mo	- seismic moment, in Nm
%	      Mw	- moment magnitude
%
% mmai, 02/08/98
% --------------

if nargin == 1
  dim = size(in); 
end
W = dim(1); L = dim(2);

if nargin<3
mu = 3.3*1e10;			% shear modulus = rigidity [N/m^2]
end
s = mean(in(:));		% average slip over fault surface [m]
A = L*W;			% fault surface [m]

Mo = mu*s*A;			% moment [Nm]
Mw = (2/3)*(log10(Mo) - 9.05);		% moment magnitude
