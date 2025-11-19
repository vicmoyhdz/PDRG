function y = rank_transform_to_Levy(x,alpha,beta,gamma,delta)

% By Victor Hernández (victorh@hi.is). July 2025    
% Rank-transform vector x to match a Levy (stable) PDF
% Used to modify the rupture velocity distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=numel(x);
% Sort input and keep index
[~, sort_idx] = sort(x(:));

% Target stable distribution parameters
% alpha = 1.4;   % Stability (0 < alpha <= 2) 1.14
% beta = -0.9;    % Skewness (-1 <= beta <= 1) -0.99
% gamma = 0.015;     % Scale (gamma > 0) 0.035
% delta = 0.82;     % Location 0.84
% Create a stable distribution object
pd = makedist('Stable', 'Alpha', alpha, 'Beta', beta, 'Gam', gamma, 'Delta', delta);

% Generate sample from stable distribution
stable_sample = random(pd, 1, N);

% Sort the stable-distributed sample
stable_sorted = sort(stable_sample);

% Map sorted stable values to the rank order of A
B = zeros(N,1);
B(sort_idx) = stable_sorted;

y=reshape(B,size(x));
end