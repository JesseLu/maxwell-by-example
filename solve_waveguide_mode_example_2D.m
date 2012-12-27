%% solve_waveguide_mode_example_2D

    omega = 0.12;
    dims = [80 80 1];
    epsilon_wg = 13-1i;
    dir = 'y+';
    mode_num = 4;
    wg_dims = [1e9 16];
    
    if dir(1) == 'y'
        wg_dims = fliplr(wg_dims);
    end

    % A rectangle that covers the entire grid.
    my_rectangle = struct('type', 'rectangle', ...
                     'position', [0 0], ...
                     'size', [1e9 1e9], ...
                     'permittivity', 1);

    % Waveguide running in the x-direction.
    my_waveguide = struct('type', 'rectangle', ...
                     'position', [0 0], ...
                     'size', wg_dims, ...
                     'permittivity', epsilon_wg);

    mu = {ones(dims), ones(dims), ones(dims)};

    epsilon = {ones(dims), ones(dims), ones(dims)};
    epsilon = add_planar(epsilon, 10, dims(3)/2, {my_rectangle, my_waveguide});

    [s_prim, s_dual] = stretched_coordinates(omega, dims, [10 10 1]);

    [beta, E, H, J] = solve_waveguide_mode(omega, s_prim, s_dual, epsilon, ...
                {[1 dims(2)/2 1], [dims(1) dims(2)/2 dims(3)]}, dir, mode_num);

    % Get ingredient matrices and vectors.
    [A1, A2, m, e, b] = maxwell_matrices(omega, s_prim, s_dual, mu, epsilon, J); 

    % Form full matrix.
    n = prod(dims);
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    A = A1 * my_diag(m.^-1) * A2 - omega^2 * my_diag(e);

    % Solve
    x = A \ b;

    % Reshape solution and plot it.
    for k = 1 : 3
        E{k} = reshape(x((k-1)*n+1 : k*n), dims);

        subplot(1, 3, k)
        imagesc(abs(E{k})'); axis equal tight;
        set(gca, 'YDir', 'normal');
        colormap jet 
    end
    snapnow;



    % solve_waveguide_mode(omega, s_prim(2:3), s_dual(2:3), eps_wg, 2);

