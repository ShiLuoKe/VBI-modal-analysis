%% f01_loadInputParams.m
% =========================================================================
% Description : Define and initialize all fundamental parameters for the
%               periodic multi-vehicle VBI scenario.
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-18
% =========================================================================

function sysParams = f01_loadInputParams()

    %% Bridge parameters  (sysParams.bridge)
    sysParams.bridge.L   = 30;
    sysParams.bridge.n   = 10;
    sysParams.bridge.mb  = 2500 * 5.3;
    sysParams.bridge.Eb  = 27.5e9;
    sysParams.bridge.Ib  = 4.98;
    sysParams.bridge.xib = 0.02;

    %% Vehicle parameters  (sysParams.vehicle)
    sysParams.vehicle.mv  = 1e5;
    sysParams.vehicle.kv  = 5e10;
    sysParams.vehicle.v   = 100;
    sysParams.vehicle.xiv = 0.2;


    %% Multi-vehicle parameters  (sysParams.fleet)
    % d_spacing - Distance between successive vehicles [m]
    % T_total   - Total simulation duration [s]
    % x0_lead   - Initial position of the leading (first) vehicle [m]
    %             x0_lead = 0  : first vehicle starts at left end of bridge
    %             x0_lead = -L : first vehicle starts one bridge-length before
    % N_max     - Maximum number of vehicles on the bridge simultaneously
    %             (computed automatically)
    sysParams.fleet.d_spacing = 3;  
    sysParams.fleet.T_total   = 10;  
    sysParams.fleet.x0_lead   = 30;    

    % Compute maximum number of vehicles on the bridge at any time
    L = sysParams.bridge.L;
    d = sysParams.fleet.d_spacing;
    sysParams.fleet.N_max = floor(L / d) + 1;

    fprintf('    Vehicle spacing   d = %.2f m\n', d);
    fprintf('    Leading vehicle start x0 = %.2f m\n', sysParams.fleet.x0_lead);
    fprintf('    Bridge span       L = %.2f m\n', L);
    fprintf('    Max vehicles on bridge N_max = %d\n', sysParams.fleet.N_max);

    %% Solver settings  (sysParams.solver)
    sysParams.solver.dt    = 0.001;
    sysParams.solver.Fs    = 1 / sysParams.solver.dt;

    % -----------------------------------------------------------------
    % (*) KEY PARAMETER : gamma
    % -----------------------------------------------------------------
    sysParams.solver.gamma = 1;  % <-- CHANGE HERE (^_^)
    fprintf('  Configuration parameter : gamma = %g\n', sysParams.solver.gamma);

end