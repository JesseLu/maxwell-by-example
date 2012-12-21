%% solve_waveguide_mode
% Find the mode of a waveguide, as well as the current excitation for it.

%% Description
% Given a 2D (or even 1D) description of the waveguiding device, 
% a corresponding waveguide mode of a given order and frequency is found.

function [] = solve_waveguide_mode(omega, s_prim, s_dual, mu, epsilon, varargin)

%% Input parameters

%% Output parameters

    %% Build operator.
    % Taken from Section 2 of Veronis, Fan, J. of Lightwave Tech., vol. 25, no. 9, Sept 2007.
%    A = -omega^2 * eps_yx + eps_yx * [Dfy; -Dfx] * inv_eps_z * [-Dby, Dbx] - ...

    %% Estimate the largest eigenvalue
    % Use only the real part of the operator.

    %% Solve for the desired mode
    % Use only the real part of the operator.

    %% Solve for the mode with the full operator

    %% Back out the current excitation


