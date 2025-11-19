function [s] = slip_vec_from_normal(N,rakes_deg)

% By Victor Hernández (victorh@hi.is). July 2025
% Calculates the slip vector given the local normal vector and rake
%
% Inputs
% N: Mx3 matrix of normal vectors (XYZ)
% rakes_deg: Mx1 vector of rake angles (degrees)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = size(N, 1);

% Default reference vector Z (0,0,1) for all rows
Z = repmat([0, 0, 1], M, 1);

% Detect nearly vertical normals to use Y instead of Z
vertical_threshold = 0.99;
dot_products = sum(N .* Z, 2);  % Mx1
is_vertical = abs(dot_products) > vertical_threshold;

% Replace Z with Y (0,1,0) where needed
Z(is_vertical, :) = repmat([0, 1, 0], sum(is_vertical), 1);

% Compute strike direction: cross(Z, N)
strike = cross(Z, N, 2);
strike = strike ./ vecnorm(strike, 2, 2);  % normalize rows

% Dip direction: cross(N, strike)
dip = cross(N, strike, 2);
dip = dip ./ vecnorm(dip, 2, 2);  % normalize rows

% Convert rake to radians
rake_rad = deg2rad(rakes_deg);

% Compute slip vectors: slip = cos(rake)*strike + sin(rake)*dip
cos_rake = cos(rake_rad);
sin_rake = sin(rake_rad);

% Element-wise multiplication and sum
s = strike .* cos_rake + dip .* sin_rake;  % Mx3

end