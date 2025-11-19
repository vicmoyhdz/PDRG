function T = fast_marching_2d(velocity, source, dx, dy)

% By Victor Hernández (victorh@hi.is). July 2025
%
% FAST_MARCHING_2D Computes travel times using the Fast Marching Method
% with anisotropic grid spacing and 8-connectivity.
%
% Inputs:
%   velocity - 2D matrix of velocities (m/s)
%   source   - [row, col] indices of the source point
%   dx, dy   - grid spacing in x (columns) and y (rows) directions
%
% Output:
%   T        - 2D matrix of travel times (s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Convert velocity to slowness
    slowness = 1 ./ velocity;

    % Initialize travel time matrix with infinities
    [nRows, nCols] = size(velocity);
    T = inf(nRows, nCols);

    % Initialize status matrix: 0 = far, 1 = narrow band, 2 = accepted
    status = zeros(nRows, nCols);

    % Set source point
    i = source(1);
    j = source(2);
    T(i, j) = 0;
    status(i, j) = 2; % accepted

    % Define neighbor offsets (8-connectivity)
    neighbors = [ ...
        -1, 0; 1, 0; 0, -1; 0, 1;  % axial
        -1,-1; -1,1; 1,-1; 1,1];   % diagonal

    % Corresponding step lengths
    step = [dy; dy; dx; dx; ...
            sqrt(dx^2 + dy^2); sqrt(dx^2 + dy^2); ...
            sqrt(dx^2 + dy^2); sqrt(dx^2 + dy^2)];

    % Initialize narrow band with neighbors of the source
    for k = 1:8
        ni = i + neighbors(k, 1);
        nj = j + neighbors(k, 2);
        if ni >= 1 && ni <= nRows && nj >= 1 && nj <= nCols
            if status(ni, nj) == 0
                T(ni, nj) = step(k) * slowness(ni, nj);
                status(ni, nj) = 1; % narrow band
            end
        end
    end

    % Main loop
    while true
        % Find narrow band points
        [nb_i, nb_j] = find(status == 1);
        if isempty(nb_i)
            break;
        end

        % Find the narrow band point with the smallest T
        nb_indices = sub2ind(size(T), nb_i, nb_j);
        [~, idx] = min(T(nb_indices));
        i = nb_i(idx);
        j = nb_j(idx);
        status(i, j) = 2; % accepted

        % Update neighbors
        for k = 1:8
            ni = i + neighbors(k, 1);
            nj = j + neighbors(k, 2);
            if ni >= 1 && ni <= nRows && nj >= 1 && nj <= nCols
                if status(ni, nj) ~= 2
                    T_new = compute_T(T, status, slowness, ni, nj, dx, dy);
                    if status(ni, nj) == 0
                        T(ni, nj) = T_new;
                        status(ni, nj) = 1; % narrow band
                    elseif T_new < T(ni, nj)
                        T(ni, nj) = T_new;
                    end
                end
            end
        end
    end
end

function T_val = compute_T(T, status, slowness, i, j, dx, dy)
    % Compute updated travel time at (i, j) using anisotropic spacing and 8-connectivity
    s = slowness(i, j);

    % Initialize large values
    Tx = inf; 
    Ty = inf;

    % Axis-aligned neighbors (for quadratic update)
    if j-1 >= 1 && status(i, j-1) == 2
        Tx = min(Tx, T(i, j-1));
    end
    if j+1 <= size(T,2) && status(i, j+1) == 2
        Tx = min(Tx, T(i, j+1));
    end
    if i-1 >= 1 && status(i-1, j) == 2
        Ty = min(Ty, T(i-1, j));
    end
    if i+1 <= size(T,1) && status(i+1, j) == 2
        Ty = min(Ty, T(i+1, j));
    end

    % Start with diagonal contributions (first-order)
    T_val = inf;
    for di = -1:2:1
        for dj = -1:2:1
            ni = i + di;
            nj = j + dj;
            if ni >= 1 && ni <= size(T,1) && nj >= 1 && nj <= size(T,2)
                if status(ni, nj) == 2
                    T_val = min(T_val, T(ni, nj) + sqrt(dx^2 + dy^2) * s);
                end
            end
        end
    end

    % Case 1: both axis neighbors are inf
    if isinf(Tx) && isinf(Ty)
        return; % diagonal update only
    end

    % Case 2: only one finite axis neighbor
    if isinf(Tx)
        T_val = min(T_val, Ty + dy * s);
        return;
    elseif isinf(Ty)
        T_val = min(T_val, Tx + dx * s);
        return;
    end

    % Case 3: both finite axis neighbors → solve quadratic
    a = 1/dx^2 + 1/dy^2;
    b = -2*Tx/dx^2 - 2*Ty/dy^2;
    c = (Tx^2)/dx^2 + (Ty^2)/dy^2 - s^2;

    disc = b^2 - 4*a*c;
    if disc < 0
        T_candidate = min(Tx + dx*s, Ty + dy*s);
    else
        T_candidate = (-b + sqrt(disc)) / (2*a);
        if T_candidate < max(Tx, Ty)
            T_candidate = min(Tx + dx*s, Ty + dy*s);
        end
    end

    T_val = min(T_val, T_candidate);
end
