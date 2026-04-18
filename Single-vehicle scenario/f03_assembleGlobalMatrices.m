%% f03_assembleGlobalMatrices.m
% =========================================================================
% Description : Assemble the time-varying global matrices (M, C, K, P)
%               for the VBI (Vehicle-Bridge Interaction) system.
%               Each time step has its own matrix, stored as a 3D array.
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-08
% -------------------------------------------------------------------------
% Input  : sysParams  - System parameter structure (from subfunction f01 & f02)
%          beamElem   - Beam element matrices (from subfunction f02)
% Output : globalMat  - Global matrices structure
%                       globalMat.M  (Dof x Dof x Num_steps)
%                       globalMat.C  (Dof x Dof x Num_steps)
%                       globalMat.K  (Dof x Dof x Num_steps)
%                       globalMat.P  (Dof x 1   x Num_steps)
% =========================================================================

function globalMat = f03_assembleGlobalMatrices(sysParams, beamElem)

    %% Unpack parameters
    L   = sysParams.bridge.L;
    n   = sysParams.bridge.n;
    v   = sysParams.vehicle.v;
    mv  = sysParams.vehicle.mv;
    cv  = sysParams.vehicle.cv;
    kv  = sysParams.vehicle.kv;
    dt  = sysParams.solver.dt;

    %% Derived quantities
    Le        = L / n;
    Dof       = 2*(n+1) + 1;           % Bridge DOFs + 1 vehicle DOF
    Dof_car   = 2*(n+1) + 1;           % Global DOF index of the vehicle
    Num_steps = round(L/v / dt) + 1;   % Number of time steps (vehicle crossing)

    %% Initialize global matrices
    globalMat.M = zeros(Dof, Dof, Num_steps);
    globalMat.C = zeros(Dof, Dof, Num_steps);
    globalMat.K = zeros(Dof, Dof, Num_steps);
    globalMat.P = zeros(Dof, 1,   Num_steps);

    fprintf('-----\n');
    fprintf('  Assembling global matrices, please wait ... (^_^)       \n');

    %% Assemble matrices for each time step
    for T = 1:Num_steps

        %% Assemble bridge element matrices
        for i = 1:n
            j   = i + 1;
            idx = 2*i-1 : 2*j;
            globalMat.M(idx, idx, T) = globalMat.M(idx, idx, T) + beamElem.Me;
            globalMat.C(idx, idx, T) = globalMat.C(idx, idx, T) + beamElem.Ce;
            globalMat.K(idx, idx, T) = globalMat.K(idx, idx, T) + beamElem.Ke;
        end

        %% Compute vehicle position at current time step
        t  = (T-1) * dt;
        xv = v * t;   % Vehicle position along the bridge [m]

        %% Assemble VBI element contribution
        if xv >= 0 && xv <= L

            % Locate element and local coordinate
            e      = min(n, floor(xv/Le) + 1);   % Element index (1..n)
            xi_loc = xv - (e-1) * Le;             % Local coordinate within element [m]

            % Hermite shape functions
            H = f03a_hermiteShape(xi_loc, Le);

            % VBI element matrices
            VBI  = f03b_VBIElement(H, sysParams);
            Mvbi = VBI.M;
            Cvbi = VBI.C;
            Kvbi = VBI.K;
            Pvbi = VBI.P;

            % Global DOF indices of the beam element currently in contact
            % with the vehicle (element e spans DOFs: 2e-1, 2e, 2e+1, 2e+2)
            gid = [2*e-1, 2*e, 2*(e+1)-1, 2*(e+1)];

            % Assemble VBI contribution into global matrices
            % --- Mass ---
            globalMat.M(Dof_car, Dof_car, T) = globalMat.M(Dof_car, Dof_car, T) + Mvbi(1,1);
            globalMat.M(Dof_car, gid,     T) = globalMat.M(Dof_car, gid,     T) + Mvbi(1, 2:5);
            globalMat.M(gid,     Dof_car, T) = globalMat.M(gid,     Dof_car, T) + Mvbi(2:5, 1);
            globalMat.M(gid,     gid,     T) = globalMat.M(gid,     gid,     T) + Mvbi(2:5, 2:5);

            % --- Damping ---
            globalMat.C(Dof_car, Dof_car, T) = globalMat.C(Dof_car, Dof_car, T) + Cvbi(1,1);
            globalMat.C(Dof_car, gid,     T) = globalMat.C(Dof_car, gid,     T) + Cvbi(1, 2:5);
            globalMat.C(gid,     Dof_car, T) = globalMat.C(gid,     Dof_car, T) + Cvbi(2:5, 1);
            globalMat.C(gid,     gid,     T) = globalMat.C(gid,     gid,     T) + Cvbi(2:5, 2:5);

            % --- Stiffness ---
            globalMat.K(Dof_car, Dof_car, T) = globalMat.K(Dof_car, Dof_car, T) + Kvbi(1,1);
            globalMat.K(Dof_car, gid,     T) = globalMat.K(Dof_car, gid,     T) + Kvbi(1, 2:5);
            globalMat.K(gid,     Dof_car, T) = globalMat.K(gid,     Dof_car, T) + Kvbi(2:5, 1);
            globalMat.K(gid,     gid,     T) = globalMat.K(gid,     gid,     T) + Kvbi(2:5, 2:5);

            % --- Load ---
            globalMat.P(gid, 1, T) = globalMat.P(gid, 1, T) + Pvbi(2:5, 1);

        else
            % Vehicle is off the bridge: add uncoupled vehicle DOF only
            globalMat.M(Dof_car, Dof_car, T) = globalMat.M(Dof_car, Dof_car, T) + mv;
            globalMat.C(Dof_car, Dof_car, T) = globalMat.C(Dof_car, Dof_car, T) + cv;
            globalMat.K(Dof_car, Dof_car, T) = globalMat.K(Dof_car, Dof_car, T) + kv;
        end

    end

    fprintf('-----\n');
    fprintf('  Global matrices assembled successfully.                  \n');
    fprintf('  Total time steps : %d\n', Num_steps);

end