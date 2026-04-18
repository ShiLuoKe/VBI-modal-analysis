%% f03b_VBIElement.m
% =========================================================================
% Description : Generate the VBI (Vehicle-Bridge Interaction) element
%               matrices for a single sprung mass on a Bernoulli-Euler beam.
%               Based on the "Present" formulation (gamma varies).
%
%               The VBI element matrices are formulated as:
%
%               [M] = [ mv          gamma*mv*{N}^T ]
%                     [ 0               [Mb]       ]
%
%               [C] = [ cv       cv*(gamma-1)*{N}^T + 2*gamma*mv*V*{N'}^T ]
%                     [-cv*{N}e  [Cb] + cv*(1-gamma)*{N}e*{N}^T           ]
%
%               [K] = [ kv       kv*(gamma-1)*{N}^T + cv*V*(gamma-1)*{N'}^T + gamma*mv*V^2*{N''}^T ]
%                     [-kv*{N}  [Kb] + kv*(1-gamma)*{N}*{N}^T + cv*V*(1-gamma)*{N}*{N'}^T      ]
%
% Author      : SLK
% Date        : 2026-04-08
% -------------------------------------------------------------------------
% Input  : H          - Shape function matrix (4 x 3) (from subfunction f03a)
%                       H(:,1) : N     Shape functions
%                       H(:,2) : dN    First derivatives
%                       H(:,3) : d2N   Second derivatives
%          sysParams  - System parameter structure (from subfunction f01 & f02)
% Output : VBI        - Structure containing VBI element matrices
%                       VBI.M  (5 x 5) Mass matrix
%                       VBI.C  (5 x 5) Damping matrix
%                       VBI.K  (5 x 5) Stiffness matrix
%                       VBI.P  (5 x 1) Load vector
% =========================================================================

function VBI = f03b_VBIElement(H, sysParams)

    %% Unpack parameters
    mv    = sysParams.vehicle.mv;
    cv    = sysParams.vehicle.cv;
    kv    = sysParams.vehicle.kv;
    v     = sysParams.vehicle.v;
    gamma = sysParams.solver.gamma;

    %% Unpack shape functions
    N   = H(:,1);   % (4 x 1)
    dN  = H(:,2);   % (4 x 1)
    d2N = H(:,3);   % (4 x 1)

    %% Initialize VBI element matrices
    VBI.M = zeros(5,5);
    VBI.C = zeros(5,5);
    VBI.K = zeros(5,5);
    VBI.P = zeros(5,1);

    %% Mass matrix  (A.1)
    VBI.M(1,1)   = mv;
    VBI.M(1,2:5) = gamma * mv * N';

    %% Damping matrix  (A.2)
    VBI.C(1,1)   = cv;
    VBI.C(1,2:5) = cv*(gamma-1)*N' + gamma*mv*2*v*dN';
    VBI.C(2:5,1) = -cv * N;
    VBI.C(2:5,2:5) = cv*(1-gamma) * (N*N');

    %% Stiffness matrix  (A.3)
    VBI.K(1,1)   = kv;
    VBI.K(1,2:5) = kv*(gamma-1)*N' + cv*v*(gamma-1)*dN' + gamma*mv*v^2*d2N';
    VBI.K(2:5,1) = -kv * N;
    VBI.K(2:5,2:5) = kv*(1-gamma)*(N*N') + cv*v*(1-gamma)*(N*dN');

    %% Load vector
    VBI.P(2:5,1) = -mv * 9.81 * N;

end