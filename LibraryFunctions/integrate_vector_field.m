function ISD = integrate_vector_field(Sx, Sy, sampling)

% By Victor Hernández (victorh@hi.is). July 2025
% This function solves the Poisson equation, i.e., it gives the scalar field (ISD) whose
% gradient best matches the Stress Drop vector field S=(Sx, Sy), where Sx and Sy
% are the stress drops along the fault's strike and dip directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pad=ceil(size(Sx,1)*0.1);
dx=sampling(2);
dy=sampling(1);
    % Pad fields
    Sx = padarray(Sx, [pad pad], 'replicate', 'both');
    Sy = padarray(Sy, [pad pad], 'replicate', 'both');

    [ny, nx] = size(Sx);

    % Compute divergence (centered differences)
    divS = zeros(ny, nx);
    divS(2:end-1, 2:end-1) = ...
        (Sx(2:end-1, 3:end) - Sx(2:end-1, 1:end-2)) / (2*dx) + ...
        (Sy(3:end, 2:end-1) - Sy(1:end-2, 2:end-1)) / (2*dy);

    % Flatten divergence
    b = reshape(divS, ny*nx, 1);

    % Build Laplacian matrix explicitly
    A = build_laplacian_matrix(nx, ny, dx, dy);

    % Fix Neumann singularity by fixing phi(1,1) = 0
    A(1,:) = 0;
    A(1,1) = 1;
    b(1) = 0;

    % Solve system
    phi_vec = A \ b;
    phi_full = reshape(phi_vec, ny, nx);

    % Remove padding
    ISD = phi_full(1+pad:end-pad, 1+pad:end-pad);
end

function A = build_laplacian_matrix(nx, ny, dx, dy)
    N = nx * ny;
    e = ones(N,1);

    % Initialize sparse matrix data holders
    rows = [];
    cols = [];
    vals = [];

    % Helper to convert 2D index (i,j) to 1D linear index
    idx = @(i,j) i + (j-1)*ny;

    for j = 1:nx
        for i = 1:ny
            k = idx(i,j);

            % Center point
            rows(end+1) = k;
            cols(end+1) = k;
            vals(end+1) = -2/dx^2 - 2/dy^2;

            % Left neighbor (if any)
            if j > 1
                rows(end+1) = k;
                cols(end+1) = idx(i,j-1);
                vals(end+1) = 1/dx^2;
            end

            % Right neighbor
            if j < nx
                rows(end+1) = k;
                cols(end+1) = idx(i,j+1);
                vals(end+1) = 1/dx^2;
            end

            % Down neighbor
            if i > 1
                rows(end+1) = k;
                cols(end+1) = idx(i-1,j);
                vals(end+1) = 1/dy^2;
            end

            % Up neighbor
            if i < ny
                rows(end+1) = k;
                cols(end+1) = idx(i+1,j);
                vals(end+1) = 1/dy^2;
            end
        end
    end

    % Build sparse matrix
    A = sparse(rows, cols, vals, N, N);
end
