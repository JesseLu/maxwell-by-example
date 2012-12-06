
%% mode_optimization_example
% Example of the optimization of an eigenmode.

%% Description
% This script varies the horizontal position of the holes of a beam resonator
% in order to match a particular resonance frequency, and to increase quality 
% factor.

% Make this a function instead of a script to allow for nested function definitions.
function [] = optimize_2D_mode_example()

%% Create the initial structure 
% We use the |add_planar| and |stretched_coordinates| functions to create our 
% structure as well as our simulation grid.
    
    spec(1) = struct(   'omega', 0.154, ...
                        'target_omega', 0.12, ...
                        'target_kappa', 0, ...
                        'polarization', 2);
    
    spec(2) = struct(   'omega', 0.254, ...
                        'target_omega', 0.24, ...
                        'target_kappa', 0, ...
                        'polarization', 3);

    dims = [200 40 1]; % Size of the simulation.
    make_structure = @(p) my_structure(dims, p);

    lattice_spacing = 12;
    p = lattice_spacing * [0.75:1:6]'; % Starting structure parameters.
    epsilon_init = make_structure(p);

    for k = 1 : length(spec)
        modes(k) = my_mode( dims, ...
                            spec(k).omega, ...
                            spec(k).target_omega, ...
                            spec(k).target_kappa, ...
                            spec(k).polarization, ...
                            make_structure, ...
                            epsilon_init); 
    end

    optimize_2D_modes(modes, p, dims, @(x) false, 15, @my_simulate, @(p, v) vis_progress(dims, [spec.polarization], p, v));
end

%% Source code for private functions

function [mode] = my_mode(dims, omega, omega_target, imag_omega_target, pol, make_structure, epsilon)


    [s_prim, s_dual] = stretched_coordinates(omega, dims, [10 10 0]); % s-parameters.

    J = {zeros(dims), zeros(dims), zeros(dims)};
    J{pol}(dims(1)/2 + 1, dims(2)/2, 1) = 1; % Central current source,

    mu = {ones(dims), ones(dims), ones(dims)}; % Permeability.

    v_guess = my_simulate(omega, s_prim, s_dual, mu, epsilon, J); % Simulate to get the guess.

    mode = struct(  'tr', omega_target, ...
                    'ti', imag_omega_target, ...
                    'v_init', v_guess, ...
                    's_prim', {s_prim}, ...
                    's_dual', {s_dual}, ...
                    'mu', {mu}, ...
                    'eig_vis', @(lambda, v) eig_vis(dims, pol, lambda, v), ...
                    'make_structure', make_structure);
end

function eig_vis(dims, pol, lambda, v)
    subplot 211; 
    n = prod(dims);
    unvec = @(z) {reshape(z(1:n), dims), reshape(z(n+1:2*n), dims), reshape(z(2*n+1:3*n), dims)};
    F = unvec(v);
    imagesc(abs(F{pol})'); axis equal tight;
    subplot 212;
end

function vis_progress(dims, pol, p, v)
    N = length(v);
    n = prod(dims);
    unvec = @(z) {reshape(z(1:n), dims), reshape(z(n+1:2*n), dims), reshape(z(2*n+1:3*n), dims)};
    function my_plot(img_data, ind) 
        subplot(N+1, 1, ind);
        imagesc(abs(img_data)'); axis equal tight;
    end
    epsilon = my_structure(dims, p);
    my_plot(epsilon{3}, 1);
    for k = 1 : N
        E = unvec(v{k});
        my_plot(E{pol(k)}, k+1);
    end 
end
 
function [epsilon] = my_structure(dims, hole_y_pos)
% Private function to create a photonic crystal beam structure.

    my_shapes = {struct('type', 'rectangle', ...
                        'position', [0 0], ...
                        'size', [1e9 1e9], ...
                        'permittivity', 1), ...
                struct('type', 'rectangle', ...
                        'position', [0 0], ...
                        'size', [1e9 12], ...
                        'permittivity', 12.25)};

    hole_radius = 4;

    for k = 1 : length(hole_y_pos)
        my_shapes{end+1} = struct('type', 'circle', ...
                                'position', [hole_y_pos(k) 0], ...
                                'radius', hole_radius, ...
                                'permittivity', 1);
        my_shapes{end+1} = struct('type', 'circle', ...
                                'position', [-hole_y_pos(k) 0], ...
                                'radius', hole_radius, ...
                                'permittivity', 1);
    end

    for k = [-1, 1]
        my_shapes{end+1} = struct('type', 'rectangle', ...
                                'position', [k*dims(1)/2 0], ...
                                'size', [20 2*hole_radius+1], ...
                                'permittivity', 12.25);
    end

    epsilon = {ones(dims), ones(dims), ones(dims)};
    epsilon = add_planar(epsilon, 1e9, 1, my_shapes);

end


function [x] = my_simulate(omega, s_prim, s_dual, mu, epsilon, J)
% Private function to simulate. Used to get initial guess.

    % Get matrices.
    [A1, A2, m, e, b] = maxwell_matrices(omega, s_prim, s_dual, mu, epsilon, J); 

    % Solve.
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    x = (A1 * my_diag(m.^-1) * A2 - omega^2 * my_diag(e)) \ b;

end

function my_plotter(x, dims)
% Private function to visualize the electric field.
    xyz = 'xyz';
    colormap jet 
    n = prod(dims);
    for k = 1 : 3
        E{k} = reshape(x((k-1)*n+1 : k*n), dims);

        subplot(3, 2, 2*k-1)
        imagesc(abs(E{k})'); 
        axis equal tight;
        title(xyz(k));
    end
    subplot(1, 2, 2);
end



function [lambda, v, w] = my_eigensolver(sim, vis, s_prim, s_dual, mu, epsilon, v_guess)
% Private function to obtain the left- and right-eigenmode of the structure.
    
    % Get ingredient matrices and vectors.
    [A1, A2, m, e] = maxwell_matrices(0, s_prim, s_dual, mu, epsilon, epsilon); 

    dims = size(epsilon{1});
    n = prod(dims);
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    unvec = @(z) {reshape(z(1:n), dims), reshape(z(n+1:2*n), dims), reshape(z(2*n+1:3*n), dims)};

    % Form full matrix.
    % Note that we actually form the matrix for F, where F = sqrt(e) E.
    A = my_diag(e.^-0.5) * A1 * my_diag(m.^-1) * A2 * my_diag(e.^-0.5);
    v_guess = sqrt(e) .* v_guess; % Convert from F-field to E-field.

    % Compose function handles.
    mult_A = @(x) A * x;

    function [x] = solve_A_shifted(lambda, b) % This is an F-field solver.
        omega = sqrt(lambda);
        J = unvec(-i * omega * b);
        x = sqrt(e) .* sim(omega, s_prim, s_dual, mu, epsilon, J);
    end
        
    % Find the eigenmode
    [lambda, v] = eigenmode_solver(mult_A, @solve_A_shifted, vis, v_guess, 10, 1e-6);
 
    % Convert v from F-field to E-field.
    v = v ./ sqrt(e);

    % Form symmetrization matrix S to obtain right-eigenmode w.  
    [spx, spy, spz] = ndgrid(s_prim{1}, s_prim{2}, s_prim{3});
    [sdx, sdy, sdz] = ndgrid(s_dual{1}, s_dual{2}, s_dual{3});

    S = my_diag([sdx(:).*spy(:).*spz(:); ...
                spx(:).*sdy(:).*spz(:); ...
                spx(:).*spy(:).*sdz(:)]);
     
    % Obtain right eigenvector.
    w = conj(S * v);
end


