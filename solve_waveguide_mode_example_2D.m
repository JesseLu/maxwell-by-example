%% solve_waveguide_mode_example_2D
function [pow, beta] = solve_waveguide_mode_example_2D(omega)
    % omega = 0.3;
    dims = [80 80 1];
    epsilon_wg = 8;
    dir = 'y+';
    mode_num = 3;
    wg_dims = [1e9 8];
    
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
    epsilon = add_planar(epsilon, 6, dims(3)/2, {my_rectangle, my_waveguide});

    [s_prim, s_dual] = stretched_coordinates(omega, dims, [0 10 0]);

    % temp = mu;
%     for k = 1 : 3
%         mu{k} = real(epsilon{k});
%     end
    % epsilon = temp;

    figure(1);
    [beta, E, H, J] = solve_waveguide_mode( ...
                omega, s_prim, s_dual, mu, epsilon, ...
                {[1 dims(2)/2 1], [dims(1) dims(2)/2 dims(3)]}, dir, mode_num);

    % Get ingredient matrices and vectors.
    [A1, A2, m, e, b] = maxwell_matrices(omega, s_prim, s_dual, mu, epsilon, J); 

    % Form full matrix.
    n = prod(dims);
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    A = A1 * my_diag(m.^-1) * A2 - omega^2 * my_diag(e);

    % Solve
    x = A \ b;

    % Get H-field.
    y = my_diag(1./(-i*m*omega)) * (A2 * x);

    % Reshape solution and plot it.
    figure(2);
    for k = 1 : 3
        E{k} = reshape(x((k-1)*n+1 : k*n), dims);
        H{k} = reshape(y((k-1)*n+1 : k*n), dims);

        subplot(2, 3, k); imagesc(real(E{k})'); axis equal tight;
        subplot(2, 3, k+3); imagesc(real(H{k})'); axis equal tight;
        set(gca, 'YDir', 'normal');
        colormap jet 
    end
    snapnow;

    % Calculate the radiant power at every y-plane.
    for k = 1 : dims(2)
        avg_Hx = 0.5 * (H{1}(:,k,1) + H{1}(:,mod(k-2,dims(2))+1,1));
        avg_Hz = 0.5 * (H{3}(:,k,1) + H{3}(:,mod(k-2,dims(2))+1,1));
        p(k) = dot(E{3}(:,k,1), avg_Hx) + dot(E{1}(:,k,1), avg_Hz);
    end
    subplot(2, 3, 4:5); plot(abs(p));

    pow = [mean(abs(p(10:30))); mean(abs(p(50:70)))];


    % solve_waveguide_mode(omega, s_prim(2:3), s_dual(2:3), eps_wg, 2);

