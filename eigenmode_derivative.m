function [df_dp] = eigenmode_derivative(lambda, v, w, make_structure, p, df_dl)

    vec = @(z) [z{1}(:); z{2}(:); z{3}(:)];

    % Compute the algebraic derivative.
    e = vec(make_structure(p));
    dl_de = -(lambda / (w' * (e .* v))) * (w' .* v.');

    % Compute the structural derivative.
    for k = 1 : length(p)
        dp = zeros(size(p));
        dp(k) = 1;
        de = 1e6 * (vec(make_structure(p + 1e-6*dp)) - e);
        dl_dp(k) = dl_de * de;
    end
    
    % Compute the objective derivative.
    df_dp = df_dl * dl_dp;

%         % Check the algebraic derivative.
%         fun = @(e) my_eigensolver(s_prim, s_dual, mu, unvec(e), v);
%         alg_err = test_derivative(fun, dl_de, lambda, e, 1e-2);
% 
%         % Check the structural derivative.
%         fun = @(p) my_eig(p, v);
%         struct_err = test_derivative(fun, dl_dp, lambda, p, 1e-2);
% 
%         % Check objective derivative.
%         fun1 = @(p) my_eig(p, v);
%         fun = @(p) f(fun1(p));
%         obj_err = test_derivative(fun, df_dp, f(lambda), p, 1e-2);
% 
%         fprintf('Derivative errors: %e, %e, %e\n', alg_err, struct_err, obj_err);

% 
% function [err] = test_derivative(fun, df_dz, f0, z0, step_len)
% % Check a derivative.
%     
%     % Produce a random direction.
%     dz = randn(size(z0));
%     dz = step_len * dz / norm(dz);
% 
%     % Evaluate delta in that direction empirically
%     delta_empirical = fun(z0 + dz) - f0;
%     delta_derivative = real(df_dz * dz);
% 
%     err = norm(delta_empirical - delta_derivative) / norm(delta_empirical);
% 
