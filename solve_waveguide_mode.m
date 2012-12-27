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
        sp{k} = s_prim{k}(p0(k):p1(k));
        sd{k} = s_dual{k}(p0(k):p1(k));
        eps{k} = epsilon{k}(p0(1):p1(1), p0(2):p1(2), p0(3):p1(3));
    end

    % Figure out what direction we are to propagate in.
    if all(dir(1) ~= 'xyz')
        error('The propagation direction must be either x, y, or z.');  
    end
    prop_dir = find(dir(1) == 'xyz');

    
    %% Build the operator
    % Build both real-only and full-complex versions of the operator.

    % Full complex operator.
    A = wg_operator(omega, sp, sd, eps, prop_dir, shape);
    n = size(A, 1);

    % Real-only operator.
    for k = 1 : 3
        sp_r{k} = real(sp{k});
        sd_r{k} = real(sd{k});
        eps_r{k} = real(eps{k});
    end
    A_r = wg_operator(real(omega), sp_r, sd_r, eps_r, prop_dir, shape);


    %% Solve for largest-magnitude eigenvalue of the real operator 
    % This is done in order to obtain the appropriate shift, 
    % from which we can calculate the most negative eigenvalues.

    % Use the power iteration algorithm.
    v = randn(n, 1);
    for k = 1 : 20 
        v = A_r * v;
    end
    ev_max = (v' * A_r * v) / norm(v)^2; % Rayleigh quotient.
    shift = abs(ev_max); 


    %% Solve for the desired eigenvector of the real operator
    % Taking the real operator, we a few of the most negative eigenmodes,
    % and then choose the one we are interested in.

    % Shift the matrix and find the appropriate eigenmode.
    % Find a few extra modes just to be sure we found the correct one.
    [V, D] = eigs(A_r - shift * speye(n), mode_num + 5); 
    
    gamma = diag(D);
    [temp, ind] = sort(gamma); % Sort most negative first.
    v = V(:,ind(mode_num)); % Choose the appropriate eigenmode.


    %% Solve for the eigenvector of the full operator
    % We use the selected eigenvector from the real operator as an initial
    % guess.

    % Perform Rayleigh quotient iteration to get the mode of the full operator.
    lambda = v' * A * v;
    for k = 1 : 40 
        err(k) = norm(A*v - lambda*v);
        if (err(k) < 1e-13)
            break
        end
        w = (A - lambda*speye(n)) \ v; 
        v = w / norm(w);
        lambda = v' * A * v;
    end


    %% Calculate output parameters

    % Wave-vector.
    beta = i * sqrt(lambda);
    % norm((A - (i*beta)^2 * speye(n))*v) % Can be used to calculate error.

    % H-field.
    N = prod(shape);
    hx = v(1:N,:);
    hy = v(N+1:end,:);

    % E-field.

    % J-field (current source excitation).

    % Plot hx and hy.
    subplot 121; my_plot(reshape(abs(hx), shape));
    subplot 122; my_plot(reshape(abs(hy), shape));

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
function [A, HtoE] = wg_operator(omega, s_prim, s_dual, epsilon, prop_dir, shape)
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
