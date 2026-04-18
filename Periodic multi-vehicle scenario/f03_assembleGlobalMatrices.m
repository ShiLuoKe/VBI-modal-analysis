%% f03_assembleGlobalMatrices.m
% =========================================================================
% Description : Assemble the time-varying global matrices (M, C, K, P)
%               for the periodic multi-vehicle VBI system.
%               At each time step, determine which vehicles are on the
%               bridge and assemble their VBI element contributions.
%
%               DOF arrangement:
%                 DOFs 1 .. 2*(n+1)              : Bridge DOFs
%                 DOFs 2*(n+1)+1 .. 2*(n+1)+N_max : Vehicle DOFs (slots)
%
%               Vehicle slot assignment:
%                 Vehicles are assigned to slots in a cyclic manner.
%                 When a vehicle leaves the bridge, its slot is released
%                 and reused by the next vehicle entering the bridge.
%
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-18
% -------------------------------------------------------------------------
% Input  : sysParams  - System parameter structure
%          beamElem   - Beam element matrices
% Output : globalMat  - Global matrices structure
% =========================================================================

function globalMat = f03_assembleGlobalMatrices(sysParams, beamElem)

    %% Unpack parameters
    L         = sysParams.bridge.L;
    n         = sysParams.bridge.n;
    v         = sysParams.vehicle.v;
    mv        = sysParams.vehicle.mv;
    cv        = sysParams.vehicle.cv;
    kv        = sysParams.vehicle.kv;
    dt        = sysParams.solver.dt;
    N_max     = sysParams.fleet.N_max;
    d_spacing = sysParams.fleet.d_spacing;
    T_total   = sysParams.fleet.T_total;
    x0_lead   = sysParams.fleet.x0_lead;

    %% Derived quantities
    Le        = L / n;
    Dof_bridge = 2*(n+1);
    Dof       = Dof_bridge + N_max;      % Total system DOFs
    Num_steps = round(T_total / dt) + 1;


    %% Determine total number of vehicles in the fleet
    %  Vehicle k position: x_k(t) = x0_lead - (k-1)*d_spacing + v*t
    %  A vehicle is needed if it can reach x = 0 before t = T_total,
    %  i.e. x0_lead - (k-1)*d_spacing + v*T_total >= 0
    %  =>   k <= (x0_lead + v*T_total) / d_spacing + 1
    N_total = ceil((x0_lead + v * T_total) / d_spacing) + 1;

    % Initial position of each vehicle (vehicle index k = 1, 2, ...)
    % x_k(t) = -(k-1)*d_spacing + v*t
    fprintf('-----\n');
    fprintf('  Total vehicles in fleet : %d\n', N_total);
    fprintf('  Total DOFs              : %d (bridge: %d + vehicle slots: %d)\n', ...
            Dof, Dof_bridge, N_max);
    fprintf('  Total time steps        : %d\n', Num_steps);

    %% Initialize global matrices
    globalMat.M = zeros(Dof, Dof, Num_steps);
    globalMat.C = zeros(Dof, Dof, Num_steps);
    globalMat.K = zeros(Dof, Dof, Num_steps);
    globalMat.P = zeros(Dof, 1,   Num_steps);

    fprintf('-----\n');
    fprintf('  Assembling global matrices, please wait ... (^_^)       \n');

    %% Assemble matrices for each time step
    for T = 1:Num_steps

        t = (T-1) * dt;

        %% Assemble bridge element matrices
        for i = 1:n
            j   = i + 1;
            idx = 2*i-1 : 2*j;
            globalMat.M(idx, idx, T) = globalMat.M(idx, idx, T) + beamElem.Me;
            globalMat.C(idx, idx, T) = globalMat.C(idx, idx, T) + beamElem.Ce;
            globalMat.K(idx, idx, T) = globalMat.K(idx, idx, T) + beamElem.Ke;
        end

        %% Determine which vehicles are on the bridge
        slot = 0;   % Counter for vehicle DOF slots
        for k = 1:N_total

            % Position of vehicle k at time t
            xv_k = x0_lead - (k-1) * d_spacing + v * t;

            if xv_k >= 0 && xv_k <= L
                slot = slot + 1;
                if slot > N_max
                    break;  % Safety: should not happen
                end

                % Global DOF index for this vehicle slot
                Dof_car = Dof_bridge + slot;

                % Locate element and local coordinate
                e      = min(n, floor(xv_k / Le) + 1);
                xi_loc = xv_k - (e-1) * Le;

                % Hermite shape functions
                H = f03a_hermiteShape(xi_loc, Le);

                % VBI element matrices
                VBI  = f03b_VBIElement(H, sysParams);
                Mvbi = VBI.M;
                Cvbi = VBI.C;
                Kvbi = VBI.K;
                Pvbi = VBI.P;

                % Global DOF indices of the beam element in contact
                gid = [2*e-1, 2*e, 2*(e+1)-1, 2*(e+1)];

                % Assemble VBI contribution
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
            end
        end

        %% Handle empty vehicle slots (vehicles not on the bridge)
        %  For slots that are not occupied at this time step, set them as
        %  uncoupled (identity-like) so the system remains well-conditioned.
        for s = (slot+1):N_max
            Dof_empty = Dof_bridge + s;
            globalMat.M(Dof_empty, Dof_empty, T) = mv;
            globalMat.C(Dof_empty, Dof_empty, T) = cv;
            globalMat.K(Dof_empty, Dof_empty, T) = kv;
        end

    end

    %% Store useful info
    globalMat.Num_steps  = Num_steps;
    globalMat.Dof        = Dof;
    globalMat.Dof_bridge = Dof_bridge;

    fprintf('-----\n');
    fprintf('  Global matrices assembled successfully.                  \n');

end