%% f04_solveEigenValue.m
% =========================================================================
% Description : Solve the eigenvalue problem for the coupled VBI system
%               at each time step to extract the coupled frequencies and
%               damping ratios of the vehicle and bridge.
%
%               Note: The extracted frequencies and damping ratios are
%               COUPLED system values, which differ from the uncoupled
%               theoretical values in theorFreq (from subfunction f02).
%
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-08
% -------------------------------------------------------------------------
% Input  : globalMat   - Global matrices structure (from subfunction f03)
%          sysParams   - System parameter structure (from subfunction f01 & f02)
%          theorFreq   - Theoretical uncoupled frequencies (from subfunction f02)
% Output : eigenResult - Structure containing coupled eigenvalue results
%                        eigenResult.xi_v   - Coupled vehicle damping ratio   [-]
%                        eigenResult.xi_b   - Coupled bridge damping ratio    [-]
%                        eigenResult.freq_v - Coupled vehicle frequency       [Hz]
%                        eigenResult.freq_b - Coupled bridge frequency        [Hz]
%                        eigenResult.t      - Time vector (same size as above)[s]
% =========================================================================

function eigenResult = f04_solveEigenValue(globalMat, sysParams, theorFreq)

    %% Unpack parameters
    Fs      = sysParams.solver.Fs;
    omegav  = 2*pi * theorFreq.vehicle.f;
    omegab1 = 2*pi * theorFreq.bridge.f1;
    omegab2 = 2*pi * theorFreq.bridge.f2;

    %% Validity check
    % This eigenvalue extraction routine is only valid when the vehicle
    % frequency is sufficiently separated from the 2nd bridge frequency.
    % If omegav >= 0.8 * omegab2, the eigenvalue indices (5 and 7) may no
    % longer correctly correspond to the vehicle and bridge modes, and
    % manual inspection of the eigenvalue results is required.
    if omegav >= 0.8 * omegab2
        error(['[f04] WARNING: omegav (%.4f rad/s) is not sufficiently ' ...
               'below omegab2 (%.4f rad/s).\n' ...
               'The automatic eigenvalue selection may be incorrect.\n' ...
               'Please inspect the eigenvalue results manually.'], ...
               omegav, omegab2);
    end

    %% Setup
    [~, ~, Num_steps] = size(globalMat.M);
    Cal_Step = 1;
    Step     = max(1, min(Cal_Step, Num_steps-1));

    % Index list for eigenvalue computation
    iList = 1:Step:(Num_steps-1);
    if iList(end) ~= (Num_steps-1)
        iList = [iList, (Num_steps-1)];
    end
    Num_calc = length(iList);

    %% Initialize output arrays
    eigenResult.xi_v   = zeros(1, Num_calc);
    eigenResult.xi_b   = zeros(1, Num_calc);
    eigenResult.freq_v = zeros(1, Num_calc);
    eigenResult.freq_b = zeros(1, Num_calc);

    fprintf('-----\n');
    fprintf('  Solving coupled eigenvalues for each time step ...      \n');
    tic;

    %% Eigenvalue computation loop
    for i = iList

        %% Apply boundary conditions (simply supported: pin both ends)
        % Left support (DOF 1)
        globalMat.K(1,:,i) = 0;   globalMat.K(:,1,i) = 0;
        globalMat.K(1,1,i) = 1;
        globalMat.M(1,:,i) = 0;   globalMat.M(:,1,i) = 0;
        globalMat.M(1,1,i) = 1;
        globalMat.C(1,:,i) = 0;   globalMat.C(:,1,i) = 0;
        globalMat.C(1,1,i) = 1;

        % Right support (second-to-last bridge DOF, excludes vehicle DOF)
        Dof_right = size(globalMat.M, 1) - 2;
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
        % [0   M] {x_dot}    [-M   0] {x}
        % [M   C] {x_ddot} = [ 0   K] {x_dot}
        A = [ zeros(size(Mu)),  Mu;
              Mu,               Cu];
        B = [-Mu,               zeros(size(Mu));
              zeros(size(Mu)),  Ku];

        %% Solve eigenvalue problem
        calc_id = find(iList == i);
        Lambda  = sort(eig(B, -A));
        OmegaVB = abs(Lambda);
        FreqVB  = OmegaVB / (2*pi);
        XiVB    = -real(Lambda) ./ OmegaVB;

        %% Extract coupled frequencies and damping ratios
        % The eigenvalues come in conjugate pairs, so each mode occupies
        % two consecutive indices. The first two pairs (indices 1-4) are
        % associated with the two constrained boundary DOFs and are
        % therefore meaningless. The first meaningful pair starts at
        % index 5 (odd indices only are selected here).
        %
        % Selection rule:
        %   omegav > omegab1 : bridge mode at index 5, vehicle mode at index 7
        %   omegav < omegab1 : vehicle mode at index 5, bridge mode at index 7
        if omegav > omegab1
            eigenResult.freq_b(1,calc_id) = FreqVB(5);
            eigenResult.xi_b(1,calc_id)   = XiVB(5);
            eigenResult.freq_v(1,calc_id) = FreqVB(7);
            eigenResult.xi_v(1,calc_id)   = XiVB(7);
        else
            eigenResult.freq_v(1,calc_id) = FreqVB(5);
            eigenResult.xi_v(1,calc_id)   = XiVB(5);
            eigenResult.freq_b(1,calc_id) = FreqVB(7);
            eigenResult.xi_b(1,calc_id)   = XiVB(7);
        end

    end

    %% Time vector
    eigenResult.t = (iList - 1) * (1/Fs);

    elapsed = toc;
    fprintf('-----\n');
    fprintf('  Coupled eigenvalue computation completed.               \n');
    fprintf('  Total elapsed time : %.3f seconds\n', elapsed);

end