%% f05_solveNewmarkBeta.m
% =========================================================================
% Description : Solve the dynamic response of the VBI system using the
%               Newmark-beta method (beta = 1/4, gamma_NM = 1/2),
%               i.e. the constant average acceleration method.
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-08
% -------------------------------------------------------------------------
% Input  : globalMat  - Global matrices structure (from subfunction f03)
%          sysParams  - System parameter structure (from subfunction f01 & f02)
% Output : dynResult  - Structure containing dynamic response
%                       dynResult.vehicle.d - Vehicle displacement  [m]
%                       dynResult.vehicle.v - Vehicle velocity      [m/s]
%                       dynResult.vehicle.a - Vehicle acceleration  [m/s^2]
%                       dynResult.bridge.d  - Bridge mid-span displacement  [m]
%                       dynResult.bridge.v  - Bridge mid-span velocity      [m/s]
%                       dynResult.bridge.a  - Bridge mid-span acceleration  [m/s^2]
% =========================================================================

function dynResult = f05_solveNewmarkBeta(globalMat, sysParams)

    %% Unpack parameters
    n  = sysParams.bridge.n;
    dt = sysParams.solver.dt;

    %% Setup
    Num_steps = size(globalMat.M, 3);
    Dof       = size(globalMat.M, 1);

    %% Initialize displacement, velocity, and acceleration arrays
    d = zeros(Dof, Num_steps);   % Displacement
    v = zeros(Dof, Num_steps);   % Velocity
    a = zeros(Dof, Num_steps);   % Acceleration

    fprintf('-----\n');
    fprintf('  Newmark-beta time integration started ...               \n');
    tic;

    %% Newmark-beta iteration
    for i = 1:Num_steps-1

        % Equivalent stiffness and load
        Ke = 4/dt^2 * globalMat.M(:,:,i) + ...
             2/dt   * globalMat.C(:,:,i) + ...
                      globalMat.K(:,:,i);

        pe = globalMat.P(:,1,i) + ...
             globalMat.M(:,:,i) * (4/dt^2 * d(:,i) + 4/dt * v(:,i) + a(:,i)) + ...
             globalMat.C(:,:,i) * (2/dt   * d(:,i) +         v(:,i));

        % Apply boundary conditions (simply supported)
        % Left support (DOF 1)
        Ke(1,:) = 0;   Ke(:,1) = 0;   Ke(1,1) = 1;
        pe(1)   = 0;

        % Right support (second-to-last bridge DOF, excludes vehicle DOF)
        Dof_right    = Dof - 2;
        Ke(Dof_right,:) = 0;   Ke(:,Dof_right) = 0;   Ke(Dof_right,Dof_right) = 1;
        pe(Dof_right)   = 0;

        % Solve for displacement at next time step
        d(:,i+1) = Ke \ pe;

        % Update acceleration and velocity
        a(:,i+1) = 4/dt^2 * (d(:,i+1) - d(:,i)) - 4/dt * v(:,i) - a(:,i);
        v(:,i+1) = v(:,i) + dt/2 * a(:,i) + dt/2 * a(:,i+1);

    end

    elapsed = toc;
    fprintf('-----\n');
    fprintf('  Newmark-beta integration completed.                     \n');
    fprintf('  Total elapsed time : %.3f seconds\n', elapsed);

    %% Extract vehicle response (last DOF)
    Dof_car = Dof;
    dynResult.vehicle.d = d(Dof_car, :);
    dynResult.vehicle.v = v(Dof_car, :);
    dynResult.vehicle.a = a(Dof_car, :);

    %% Extract bridge mid-span response (DOF at n + 1)
    % Valid only when n is even
    if mod(n,2) ~= 0
        warning('[f05] n is odd: no exact mid-span node exists. Using closest node.');
    end
    Dof_mid = n + 1;
    dynResult.bridge.d = d(Dof_mid, :);
    dynResult.bridge.v = v(Dof_mid, :);
    dynResult.bridge.a = a(Dof_mid, :);

end