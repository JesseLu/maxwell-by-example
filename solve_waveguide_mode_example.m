%% solve_waveguide_mode_example

    omega = 0.12;
    dims = [80 40 40];


    % A rectangle that covers the entire grid.
    my_rectangle = struct('type', 'rectangle', ...
                     'position', [0 0], ...
                     'size', [1e9 1e9], ...
                     'permittivity', 1);

    % Waveguide running in the x-direction.
    my_waveguide = struct('type', 'rectangle', ...
                     'position', [0 0], ...
                     'size', [1e9 16], ...
                     'permittivity', 13);

    epsilon = {ones(dims), ones(dims), ones(dims)};
    epsilon = add_planar(epsilon, 10, dims(3)/2, {my_rectangle, my_waveguide});

    [s_prim, s_dual] = stretched_coordinates(omega, dims, [0 0 10]);

%     % Visualize the structure.
%     for k = 1 : 3
%        subplot(2, 3, k);
%        imagesc(epsilon{k}(:,:,dims(3)/2)'); axis equal tight;
%        subplot(2, 3, k+3);
%        imagesc(squeeze(epsilon{k}(dims(1)/2,:,:))'); axis equal tight;
%     end
%     snapnow;

    for k = 1 : 3
        j = mod(k+1, 3) + 1;
        eps_wg{k} = epsilon{j}(1,:,:);
    end
    solve_waveguide_mode(omega, s_prim, s_dual, epsilon, {[20 1 1], [20 dims(2) dims(3)]}, 'x', 2)
    % solve_waveguide_mode(omega, s_prim(2:3), s_dual(2:3), eps_wg, 2);

