%% f04_solveEigenValue.m
% =========================================================================
% Description : Solve the eigenvalue problem for the coupled multi-vehicle
%               VBI system at each time step. Extract the coupled bridge
%               frequency and damping ratio.
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-18
% -------------------------------------------------------------------------
% Input  : globalMat   - Global matrices structure
%          sysParams   - System parameter structure
%          theorFreq   - Theoretical uncoupled frequencies
% Output : eigenResult - Coupled eigenvalue results (bridge mode only)
% =========================================================================

function eigenResult = f04_solveEigenValue(globalMat, sysParams, theorFreq)

    %% Unpack parameters
    Fs      = sysParams.solver.Fs;
    omegab1 = 2*pi * theorFreq.bridge.f1;

    %% Setup
    Num_steps = globalMat.Num_steps;
    Cal_Step  = 1;
    Step      = max(1, min(Cal_Step, Num_steps-1));

    iList = 1:Step:(Num_steps-1);
    if iList(end) ~= (Num_steps-1)
        iList = [iList, (Num_steps-1)];
    end
    Num_calc = length(iList);

    %% Initialize output arrays
    eigenResult.freq_b = zeros(1, Num_calc);
    eigenResult.xi_b   = zeros(1, Num_calc);

    fprintf('-----\n');
    fprintf('  Solving coupled eigenvalues for each time step ...      \n');
    tic;

    %% Eigenvalue computation loop
    for i = iList

        %% Apply boundary conditions (simply supported)
        % Left support (DOF 1)
        globalMat.K(1,:,i) = 0;   globalMat.K(:,1,i) = 0;
        globalMat.K(1,1,i) = 1;
        globalMat.M(1,:,i) = 0;   globalMat.M(:,1,i) = 0;
        globalMat.M(1,1,i) = 1;
        globalMat.C(1,:,i) = 0;   globalMat.C(:,1,i) = 0;
        globalMat.C(1,1,i) = 1;

        % Right support
        Dof_right = globalMat.Dof_bridge - 1;
        globalMat.K(Dof_right,:,i) = 0;   globalMat.K(:,Dof_right,i) = 0;
        globalMat.K(Dof_right,Dof_right,i) = 1;
        globalMat.M(Dof_right,:,i) = 0;   globalMat.M(:,Dof_right,i) = 0;
        globalMat.M(Dof_right,Dof_right,i) = 1;
        globalMat.C(Dof_right,:,i) = 0;   globalMat.C(:,Dof_right,i) = 0;
        globalMat.C(Dof_right,Dof_right,i) = 1;

        %% Extract matrices at current time step
        Mu = globalMat.M(:,:,i);
        Cu = globalMat.C(:,:,i);
        Ku = globalMat.K(:,:,i);

        %% State-space formulation
        A = [ zeros(size(Mu)),  Mu;
              Mu,               Cu];
        B = [-Mu,               zeros(size(Mu));
              zeros(size(Mu)),  Ku];

        %% Solve eigenvalue problem
        calc_id = find(iList == i);
        Lambda  = eig(B, -A);
        OmegaVB = abs(Lambda);
        FreqVB  = OmegaVB / (2*pi);
        XiVB    = -real(Lambda) ./ OmegaVB;

        %% Extract bridge first mode
        %  Strategy 1: find the eigenvalue whose frequency is closest to
        %  the theoretical uncoupled bridge frequency omegab1.
        target_fb = theorFreq.bridge.f1;
        [~, idx_b] = min(abs(FreqVB - target_fb));

        eigenResult.freq_b(1, calc_id) = FreqVB(idx_b);
        eigenResult.xi_b(1, calc_id)   = XiVB(idx_b);

%         %% Strategy 2: Manual — fixed eigenvalue index
%         eigenResult.freq_b_manual(1, calc_id) = FreqVB(manual_index);
%         eigenResult.xi_b_manual(1, calc_id)   = XiVB(manual_index);


    end

    %% Time vector
    eigenResult.t = (iList - 1) * (1/Fs);

    elapsed = toc;
    fprintf('-----\n');
    fprintf('  Coupled eigenvalue computation completed.               \n');
    fprintf('  Total elapsed time : %.3f seconds\n', elapsed);

end