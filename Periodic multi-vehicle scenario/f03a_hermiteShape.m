%% f03a_hermiteShape.m
% =========================================================================
% Description : Compute the Hermite shape functions and their derivatives
%               for a Bernoulli-Euler beam element.
%
%               DOF order: [v_i, theta_i, v_j, theta_j]
%               where v = transverse displacement, theta = rotation
%               i = left node (x=0), j = right node (x=Le)
%
% Author      : SLK
% Date        : 2026-04-08
% -------------------------------------------------------------------------
% Input  : xi_loc  - Local coordinate within the element [m]
%          Le      - Element length [m]
% Output : H       - Shape function matrix (4 x 3)
%                    H(:,1) : Shape functions        N
%                    H(:,2) : First derivatives      dN/dx
%                    H(:,3) : Second derivatives     d2N/dx2
% =========================================================================

function H = f03a_hermiteShape(xi_loc, Le)

    %% Local coordinate ratio
    s = xi_loc / Le;   % Normalised local coordinate [-]

    %% Shape functions  N
    N1 =  1 - 3*s^2 + 2*s^3;
    N2 =  xi_loc * (1 - s)^2;
    N3 =  3*s^2 - 2*s^3;
    N4 =  xi_loc^2/Le * (s - 1);
    N   = [N1; N2; N3; N4];

    %% First derivatives  dN/dx
    dN1 = -6*s/Le  + 6*s^2/Le;
    dN2 =  1 - 4*s + 3*s^2;
    dN3 =  6*s/Le  - 6*s^2/Le;
    dN4 = -2*s     + 3*s^2;
    dN  = [dN1; dN2; dN3; dN4];

    %% Second derivatives  d2N/dx2
    d2N1 = -6/Le^2 + 12*s/Le^2;
    d2N2 = -4/Le   +  6*s/Le;
    d2N3 =  6/Le^2 - 12*s/Le^2;
    d2N4 = -2/Le   +  6*s/Le;
    d2N  = [d2N1; d2N2; d2N3; d2N4];

    %% Output: (4 x 3) matrix, each column is a set of shape functions
    H = [N, dN, d2N];

end