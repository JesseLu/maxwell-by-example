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

    % Cut out relevant parameters.
    for k = 1 : 3
        s_prim{k} = s_prim{k}(p0(k):p1(k));
        s_dual{k} = s_dual{k}(p0(k):p1(k));
        epsilon{k} = epsilon{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
    end

    prop_dir
    shape = size(eps{1}); % The shape of eps 

    % Indices of the non-propagating directions.
    xdir = mod(prop_dir + 1 - 1, 3) + 1; 
    ydir = mod(prop_dir + 2 - 1, 3) + 1; 
    %% Build operators

    % Full complex operator.
    A = wg_operator(omega, s_prim, s_dual, epsilon);

    % Real-only operator.
    for k = 1 : 2
        s_prim_r{k} = real(s_prim{k});
        s_dual_r{k} = real(s_dual{k});
    end

    for k = 1 : 3
        epsilon_r{k} = real(epsilon{k});
    end
    A_r = wg_operator(real(omega), s_prim_r, s_dual_r, epsilon_r);


    %% Estimate the largest eigenvalue
    % Use only the real part of the operator.
    [V, ev_max] = eigs(A_r, 1); % Find the largest-magnitude eigenvector.

    % This shift will make the most negative eigenvectors the largest amplitude 
    % eigenvectors.
    shift = abs(ev_max); 

    %% Solve for the desired mode
    % Use only the real part of the operator.
    [V, D] = eigs(A_r - shift * speye(size(A_r,1)), mode_num); 

    %% Solve for the mode with the full operator

    gamma = diag(D);
    [temp, ind] = sort(gamma);
    beta = gamma(ind(mode_num));
    beta = 1/i * sqrt(beta + ev_max);
    V = V(:,ind(mode_num));

    % Back out hx and hy.
    shape = size(squeeze(epsilon{1}));
    size(V)
    shape
    N = prod(shape);
    hx = V(1:N,:);
    hy = V(N+1:end,:);

    % Plot hx and hy.
    for k = 1 : size(V,2)
        subplot 121; my_plot(reshape(abs(hx(:,k)), shape));
        subplot 122; my_plot(reshape(abs(hy(:,k)), shape));
    end
    beta

    % Get the current sources (notice the switch).
    jx = reshape(hy, shape);
    jy = reshape(hx, shape);


    %% Back out the current excitation



end % End of solve_waveguide_mode function.

function my_plot(x)
    imagesc(squeeze(x).', (max(abs(x(:))) + eps) * [-1 1]);
    colorbar 
    axis equal tight;
    set(gca, 'YDir', 'normal');
    end


%% Source code for private functions
function [A] = wg_operator(omega, s_prim, s_dual, epsilon)
% Taken from Section 2 of Veronis, Fan, J. of Lightwave Tech., vol. 25, no. 9, 
% Sept 2007.

    shape = size(epsilon{1});
    if length(shape) ~= 3
        error('Could not form a three-dimensional, flat shape parameter.');
    end


    % Derivative matrices.
    Dx = deriv('x', shape);
    Dy = deriv('y', shape);

    % Expand stretched-coordinate parameters.
    [s_prim_x, s_prim_y] = ndgrid(s_prim{1}, s_prim{2});
    [s_dual_x, s_dual_y] = ndgrid(s_dual{1}, s_dual{2});

    my_diag = @(z) spdiags(z(:), 0, numel(z), numel(z)); % Helper function.

    % Build the component matrices.
    Dfx = my_diag(s_dual_x) * Dx;
    Dbx = my_diag(s_prim_x) * (-Dx');
    Dfy = my_diag(s_dual_y) * Dy;
    Dby = my_diag(s_prim_y) * (-Dy');
    eps_xy = my_diag([epsilon{1}(:); epsilon{2}(:)]);
    inv_eps_z = my_diag(epsilon{3}.^-1);

    % Assemble into final operator
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
