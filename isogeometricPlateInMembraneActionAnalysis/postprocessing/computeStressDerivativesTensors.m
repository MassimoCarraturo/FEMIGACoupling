%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DsigmaDx_matrix,DsigmaDy_matrix] = get_stress_derivatives_tensors(i,p,u,U,j,q,v,V,CP,D,d)
%% Function documentation
%
% Returns the derivatives of the stress tensor with respect to the physical
% domain x-y in the parametric coordinates (u.v)
%
%                 Input : 
%                   p,q : polynomial degrees
%                   U,V : knot vectors in u,v-direction
%                    CP : Control Point coordinates and weights
%                   u,v : Parametric coordinate of the NURBS surface
%                     D : The material matrix for the given membrane
%                     d : element displacement vector
%
%                Output :
% displacement_gradient : element strain vector strain vector epsilon = [exx eyy 2exy]'
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the basis functions and its derivatives
%
% 2. Compute the Jacobian of the transformation from the physical space to the
%    NURBS parameter domain:
%
% 3. Compute the derivatives of the NURBS basis functions on the physical
%    space from those on the parent domain
%
% 4. Compute the B-operator matrices for both components of the displacement
%    fields
%
% 5. Compute the displacement gradient matix
%
%% Function main body

%% 0. Read input

% Initialize output arrays
DsigmaDx_matrix = zeros(3,3);
DsigmaDy_matrix = zeros(3,3);

% Number of basis functions per element
ne = (p+1)*(q+1);  

% Number of degrees of freedom per element
ndof = 2*ne;

% 1. Compute the basis functions and its derivatives
[~,dR,ddR] = nurbs_basis_functions_first_and_second_derivatives_2D(i,p,u,U,j,q,v,V,CP);    

%% 2. Compute the Hessian matrix which is needed for the tranformation of the second derivatives

% Initialize Jacobian matrix
Jxxi = zeros(2,2);

% Initialize Hessian matrix
Hxxi = zeros(2,3);

% initialize counter
k = 0;

% Loop over all the non-zero contributions at the span
% under study
for c = 0:q
    for b = 0:p
        % Update counter
        k = k + 1;

        % Compute recursively the entries of the Jacobian
        Jxxi(1,1) = Jxxi(1,1) + CP(i-p+b,j-q+c,1)*dR(k,1);
        Jxxi(1,2) = Jxxi(1,2) + CP(i-p+b,j-q+c,2)*dR(k,1);
        Jxxi(2,1) = Jxxi(2,1) + CP(i-p+b,j-q+c,1)*dR(k,2);
        Jxxi(2,2) = Jxxi(2,2) + CP(i-p+b,j-q+c,2)*dR(k,2);

        % Compute recursively the entries of the
        % Hessian matrix
        % d^2 x1/du^2
        Hxxi(1,1) = Hxxi(1,1) + CP(i-p+b,j-q+c,1)*ddR(k,1);
        % d^2 x1/dv^2
        Hxxi(1,2) = Hxxi(1,2) + CP(i-p+b,j-q+c,1)*ddR(k,2);
        % d^2 x1/dudv
        Hxxi(1,3) = Hxxi(1,3) + CP(i-p+b,j-q+c,1)*ddR(k,3);
        % d^2 x2/du^2
        Hxxi(2,1) = Hxxi(2,1) + CP(i-p+b,j-q+c,2)*ddR(k,1);
        % d^2 x2/dv^2
        Hxxi(2,2) = Hxxi(2,2) + CP(i-p+b,j-q+c,2)*ddR(k,2);
        % d^2 x2/dudv
        Hxxi(2,3) = Hxxi(2,3) + CP(i-p+b,j-q+c,2)*ddR(k,3);
    end
end

% On the computation of the second derivatives of the 
% basis functions in the physical space given those at
% the parent domain

% Initialize tranformation matrix for the second
% derivatives
DJxxi = zeros(5,5);

% Assign the values of DJxxi
% First block of DJxxi
for ki=1:2
    for kj=1:2
        DJxxi(ki,kj) = Jxxi(ki,kj)^2;
    end
end

% Complete first row
DJxxi(1,3) = 2*Jxxi(1,1)*Jxxi(1,2);
DJxxi(1,4) = Hxxi(1,1)^2;
DJxxi(1,5) = Hxxi(2,1)^2;

% Complete second row
DJxxi(2,3) = 2*Jxxi(2,1)*Jxxi(2,2);
DJxxi(2,4) = Hxxi(1,2)^2;
DJxxi(2,5) = Hxxi(2,2)^2;

% Complete third row
DJxxi(3,1) = Jxxi(1,1)*Jxxi(2,1);
DJxxi(3,2) = Jxxi(1,2)*Jxxi(2,2);
DJxxi(3,3) = Jxxi(1,1)*Jxxi(2,2)+Jxxi(2,1)*Jxxi(1,2);
DJxxi(3,4) = Hxxi(1,3);
DJxxi(3,5) = Hxxi(2,3);

% Last block is the Jacobian matrix
for ki=1:2
    for kj=1:2
        DJxxi(ki+3,kj+3) = Jxxi(ki,kj);
    end
end

% Initialization of the matrix with the derivatives
ddRddx = zeros(ne,3);

% loop over all the elements
for ki=1:ne
    deriv_R_actual = [ddR(ki,1:3) dR(ki,1:2)];
    deriv_Rdx = DJxxi\deriv_R_actual';
    ddRddx(ki,1:3) = deriv_Rdx(1:3);
end   

%% 3. Compute the B-operator matrices corresponding to the derivatives of the strain field with respect to the physical space

% initialize matrix
BH_1 = zeros(3,ndof);
BH_2 = zeros(3,ndof);

% Assign the values of BH_1 and BH_2 entries recursively
for i=1:ne
    % BH_1:
    % Related to d(epsilon_11)/d(dx)
    BH_1(1,2*i-1) = ddRddx(i,1);
    % Related to d(epsilon_22)/d(dx)
    BH_1(2,2*i) = ddRddx(i,3);
    % Related to d(2*epsilon_12)/d(dx)
    BH_1(3,2*i-1) = ddRddx(i,3);
    BH_1(3,2*i) = ddRddx(i,1);
    
    % BH_2:
    BH_2(1,2*i-1) = ddRddx(i,3);
    % Related to d(epsilon_22)/d(dy)
    BH_2(2,2*i) = ddRddx(i,2);
    % Related to d(2*epsilon_12)/d(dy)
    BH_2(3,2*i-1) = ddRddx(i,2);
    BH_2(3,2*i) = ddRddx(i,3);
end

%% 4. Compute the derivatives of the stress tensor using the material law

% derivative d(sigma)/(dx)
DsigmaDx = D*BH_1*d;
% derivative d(sigma)/(dx)
DsigmaDy = D*BH_2*d;

% Arrangement into matrices
% derivative d(sigma)/(dx)
DsigmaDx_matrix(1,1) = DsigmaDx(1,1);
DsigmaDx_matrix(2,2) = DsigmaDx(2,1);
DsigmaDx_matrix(1,2) = DsigmaDx(3,1);
% Due to the symmetry of the stress tensor:
DsigmaDx_matrix(2,1) = DsigmaDx(3,1);

% derivative d(sigma)/(dx)
DsigmaDy_matrix(1,1) = DsigmaDy(1,1);
DsigmaDy_matrix(2,2) = DsigmaDy(2,1);
DsigmaDy_matrix(1,2) = DsigmaDy(3,1);
% Due to the symmetry of the stress tensor:
DsigmaDy_matrix(2,1) = DsigmaDy(3,1);


end