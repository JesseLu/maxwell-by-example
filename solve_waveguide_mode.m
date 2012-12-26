%% solve_waveguide_mode
% Find the mode of a waveguide, as well as the current excitation for it.

%% Description
% Given a 2D (or even 1D) description of the waveguiding device, 
% a corresponding waveguide mode of a given order and frequency is found.

function [J, beta, E, H] = solve_waveguide_mode(...
                                        omega, s_prim, s_dual, epsilon, ...
                                        pos, dir, mode_num)

%% Input parameters
% The input parameters are very similar to those which describe a simulation,
% with the exception that most of the parameters are in two-dimensions (x and y)
% only.

%% Output parameters

    %% Parse inputs

    p0 = pos{1};
    p1 = pos{2};
    shape = p1 - p0 + 1;

    % Cut out relevant parameters.
    for k = 1 : 3
        s_prim{k} = s_prim{k}(p0(k):p1(k));
        s_dual{k} = s_dual{k}(p0(k):p1(k));
        epsilon{k} = epsilon{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
    end


    if all(dir(1) ~= 'xyz')
        error('The propagation direction must be either x, y, or z.');  
    end

    prop_dir = find(dir(1) == 'xyz');
    
    % Full complex operator.
    A = wg_operator(omega, s_prim, s_dual, epsilon, prop_dir, shape);

    % Real-only operator.
    for k = 1 : 3
        s_prim_r{k} = real(s_prim{k});
        s_dual_r{k} = real(s_dual{k});
        epsilon_r{k} = real(epsilon{k});
    end
    A_r = wg_operator(real(omega), s_prim_r, s_dual_r, epsilon_r, ...
                        prop_dir, shape);

    % Solve for largest-magnitude eigenmode.
    % This will almost always be positive.
    % We do this in order to be able to solve for the most negative eigenmodes,
    % from which we can select the appropriate propagating mode.
    [V, ev_max] = eigs(A_r, 1); % Find the largest-magnitude eigenvector.
    if ev_max < 0 
        error('Largest eigenvalue is negative.'); % Ooops.
    end

    % Shift the matrix and find the appropriate eigenmode.
    n = size(A_r, 1);
    [V, D] = eigs(A_r - ev_max * speye(n), 2*mode_num + 5); 

    gamma = diag(D);
    [temp, ind] = sort(gamma);
    beta = gamma(ind(mode_num));
    beta = 1/i * sqrt(beta + ev_max);
    V = V(:,ind(mode_num));

    % Calculate error.
    norm((A_r - beta^2 * speye(n))*V)

    % Back out hx and hy.
    N = prod(shape);
    hx = V(1:N,:);
    hy = V(N+1:end,:);

    % Plot hx and hy.
%     figure(1);
    for k = 1 : size(V,2)
        subplot 121; my_plot(reshape(real(hx(:,k)), shape));
        subplot 122; my_plot(reshape(real(hy(:,k)), shape));
    end

%     figure(2);
%     for k = 1 : 3
%         subplot(1, 3, k);
%         my_plot(epsilon{k});
%     end

    % Get the current sources (notice the switch).
    jx = reshape(hy, shape);
    jy = reshape(hx, shape);


end % End of solve_waveguide_mode function.

function my_plot(x)
    imagesc(squeeze(x).', (max(abs(x(:))) + eps) * [-1 1]);
    colorbar 
    axis equal tight;
    set(gca, 'YDir', 'normal');
    end


%% Source code for private functions
function [A] = wg_operator(omega, s_prim, s_dual, epsilon, prop_dir, shape)
% Taken from Section 2 of Veronis, Fan, J. of Lightwave Tech., vol. 25, no. 9, 
% Sept 2007.

    % Indices of the non-propagating directions.
    xdir = mod(prop_dir + 1 - 1, 3) + 1; 
    ydir = mod(prop_dir + 2 - 1, 3) + 1; 

    % Create matrices.
    xyz = 'xyz';
    Dx = deriv(xyz(xdir), shape);
    Dy = deriv(xyz(ydir), shape);

    % SC parameters.
    % Only real part is taken, this may be modified in the future.
    [s_prim_x, s_prim_y] = ndgrid(s_prim{xdir}, s_prim{ydir});
    [s_dual_x, s_dual_y] = ndgrid(s_dual{xdir}, s_dual{ydir});

    % Build matrices.
    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z));
    Dfx = my_diag(s_dual_x) * Dx;
    Dbx = my_diag(s_prim_x) * (-Dx');
    Dfy = my_diag(s_dual_y) * Dy;
    Dby = my_diag(s_prim_y) * (-Dy');
    eps_xy = my_diag([epsilon{xdir}(:); epsilon{ydir}(:)]);
    inv_eps_z = my_diag(epsilon{prop_dir}.^-1);

    % Build operator.
    A = -omega^2 * eps_xy + eps_xy * [Dfy; -Dfx] * inv_eps_z * [-Dby, Dbx] - ...
        [Dbx; Dby] * [Dfx, Dfy];

end % End of wg_operator private function.


function [D] = deriv(dir, shape)
% Private function for creating derivative matrices.
% Note that we are making the forward derivative only.
% Also, we assume periodic boundary conditions.

    shift = (dir == 'xyz'); % Direction of shift.

    % Get the displaced spatial markers.
    my_disp = @(n, shift) mod([1:n] + shift - 1, n) + 1;
    [i, j, k] = ndgrid(my_disp(shape(1), shift(1)), ...
                        my_disp(shape(2), shift(2)), ...
                        my_disp(shape(3), shift(3)));

    % Translate spatial indices into matrix indices.
    N = prod(shape);
    i_ind = 1 : N;
    j_ind = i + (j-1) * shape(1) + (k-1) * shape(1) * shape(2);

    % Create the sparse matrix.
    D = sparse([i_ind(:); i_ind(:)], [i_ind(:), j_ind(:)], ...
                [-ones(N,1); ones(N,1)], N, N);
end % End of deriv private function.
