%% f02_buildElementMatrices.m
% =========================================================================
% Description : Construct the Bernoulli-Euler beam element matrices (M, K, C)
%               and compute the theoretical uncoupled frequencies of the
%               bridge and vehicle.
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-08
% -------------------------------------------------------------------------
% Input  : sysParams  - System parameter structure (from subfunction f01)
% Output : beamElem   - Beam element matrices (Me, Ke, Ce)
%          theorFreq  - Theoretical uncoupled frequencies (bridge & vehicle)
%          sysParams  - Updated system parameter structure
%                       (adds sysParams.vehicle.cv)
% =========================================================================

function [beamElem, theorFreq, sysParams] = f02_buildElementMatrices(sysParams)

    %% Unpack parameters
    L   = sysParams.bridge.L;
    n   = sysParams.bridge.n;
    mb  = sysParams.bridge.mb;
    Eb  = sysParams.bridge.Eb;
    Ib  = sysParams.bridge.Ib;
    xib = sysParams.bridge.xib;
    mv  = sysParams.vehicle.mv;
    kv  = sysParams.vehicle.kv;
    xiv = sysParams.vehicle.xiv;
    v   = sysParams.vehicle.v;

    %% Vehicle natural frequency
    omegav = sqrt(kv / mv);
    theorFreq.vehicle.f = omegav / (2*pi);

    fprintf('-----\n');
    fprintf('  Vehicle natural frequency:\n');
    fprintf('    fv = %.4f Hz\n', theorFreq.vehicle.f);

    %% Bridge theoretical natural frequencies
    omegab1 = pi^2 / L^2 * sqrt(Eb * Ib / mb);   % 1st mode
    omegab2 = 4 * omegab1;                         % 2nd mode
    omegab3 = 9 * omegab1;                         % 3rd mode
    theorFreq.bridge.f1 = omegab1 / (2*pi);
    theorFreq.bridge.f2 = omegab2 / (2*pi);
    theorFreq.bridge.f3 = omegab3 / (2*pi);

    fprintf('-----\n');
    fprintf('  Bridge theoretical natural frequencies:\n');
    fprintf('    fb1 = %.4f Hz\n', theorFreq.bridge.f1);
    fprintf('    fb2 = %.4f Hz\n', theorFreq.bridge.f2);
    fprintf('    fb3 = %.4f Hz\n', theorFreq.bridge.f3);

    %% Beam element length
    Le = L / n;

    %% Beam element mass matrix
    beamElem.Me = mb * Le / 420 * ...
        [ 156,      22*Le,    54,     -13*Le;
          22*Le,  4*Le^2,   13*Le,  -3*Le^2;
          54,      13*Le,   156,     -22*Le;
         -13*Le, -3*Le^2,  -22*Le,   4*Le^2];

    %% Beam element stiffness matrix
    beamElem.Ke = 2*Eb*Ib / Le^3 * ...
        [  6,     3*Le,   -6,    3*Le;
           3*Le,  2*Le^2, -3*Le,  Le^2;
          -6,    -3*Le,    6,   -3*Le;
           3*Le,   Le^2,  -3*Le, 2*Le^2];

    %% Rayleigh damping coefficients
    xi  = [xib; xib];   % Assuming equal damping ratio for the first two modes
    A   = 2*omegab1*omegab2 / (omegab2^2 - omegab1^2);
    Mat = [  omegab2,   -omegab1;
            -1/omegab2,  1/omegab1 ];
    a   = A * Mat * xi;
    a0  = a(1);   % Mass-proportional Rayleigh coefficient
    a1  = a(2);   % Stiffness-proportional Rayleigh coefficient

    %% Beam element damping matrix
    beamElem.Ce = a0 * beamElem.Me + a1 * beamElem.Ke;

    %% Vehicle damping coefficient
    sysParams.vehicle.cv = 2 * mv * omegav * xiv;

    fprintf('-----\n');
    fprintf('  Vehicle damping coefficient:\n');
    fprintf('    cv = %.4f N·s/m\n', sysParams.vehicle.cv);

    %% Non-dimensional parameters
    alpha_m  = mv / (mb * L);
    alpha_f  = theorFreq.vehicle.f / theorFreq.bridge.f1;
    alpha_xi = xiv / xib;
    % Note: alpha_V is the speed parameter. In some earlier references this
    %       parameter was denoted as S. In this code and the associated
    %       paper, alpha is used consistently throughout.
    alpha_V  = pi * v / (omegab1 * L);

    fprintf('-----\n');
    fprintf('  Non-dimensional parameters:\n');
    fprintf('    Mass ratio      alpha_m  = %.4f\n', alpha_m);
    fprintf('    Frequency ratio alpha_f  = %.4f\n', alpha_f);
    fprintf('    Damping ratio   alpha_xi = %.4f\n', alpha_xi);
    fprintf('    Speed parameter alpha_V  = %.4f\n', alpha_V);

end