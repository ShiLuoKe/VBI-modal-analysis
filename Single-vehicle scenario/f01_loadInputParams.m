%% f01_loadInputParams.m
% =========================================================================
% Description : Define and initialize all fundamental parameters.
%               Controls the key variable settings used in subsequent
%               subfunctions.
% Author      : Shi LuoKe (史罗克) 
% Date        : 2026-04-08
% -------------------------------------------------------------------------
% Terminology :
%               "Bridge"  refers to the Bernoulli-Euler beam structure.
%               "Vehicle" refers to the single sprung mass moving on the
%               bridge.
% =========================================================================

function sysParams = f01_loadInputParams()

    %% Bridge parameters  (sysParams.bridge)
    % L   - Bridge span length [m]
    % n   - Number of beam elements [-]
    % mb  - Mass per unit length of the bridge [kg/m]
    % Eb  - Elastic modulus of the bridge material [Pa]
    % Ib  - Second moment of area of the bridge cross-section [m^4]
    % xib - Bridge damping ratio [-]
    sysParams.bridge.L   = 30;
    sysParams.bridge.n   = 10;
    sysParams.bridge.mb  = 2500 * 5.3;
    sysParams.bridge.Eb  = 27.5e9;
    sysParams.bridge.Ib  = 4.98;
    sysParams.bridge.xib = 0.02;

    %% Vehicle parameters  (sysParams.vehicle)
    % mv  - Vehicle mass [kg]
    % kv  - Vehicle stiffness [N/m]
    % v   - Vehicle speed [m/s]
    % xiv - Vehicle damping ratio [-]
    sysParams.vehicle.mv  = 10000;
    sysParams.vehicle.kv  = 3e7;
    sysParams.vehicle.v   = 3;
    sysParams.vehicle.xiv = 0.2;

    %% Solver settings  (sysParams.solver)
    % dt    - Integration time step [s]
    % Fs    - Sampling frequency [Hz]
    % gamma - Configuration definition parameter [-]
    %         Controls the Configuration selection of the system.
    %         Modify this value to investigate its effect on eigenvalues.
    %         gamma = 0 : Conventional coordinate selection (default)
    %         gamma = 1 : Proposed coordinate selection (recommended)
    sysParams.solver.dt    = 0.001;
    sysParams.solver.Fs    = 1 / sysParams.solver.dt;

    % -----------------------------------------------------------------
    % (*) KEY PARAMETER : Modify gamma here to switch coordinates.
    %     gamma = 0  -->  Conventional  
    %     gamma = 1  -->  Proposed      (recommended)
    % -----------------------------------------------------------------
    sysParams.solver.gamma = 0;  % <-- CHANGE HERE (^_^)
    fprintf('  Configuration parameter : gamma = %g\n', sysParams.solver.gamma);

end