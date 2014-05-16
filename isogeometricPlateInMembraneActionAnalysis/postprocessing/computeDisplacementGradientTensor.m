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
function displacement_gradient = get_displacement_gradient_tensor(i,p,u,U,j,q,v,V,CP,d)
%% Function documentation
%
% Returns the displacement gradient vector = [ux,x ux,y ; uy,x uy,y]' at
% the parametric location (u,v)
%
%                 Input : 
%                   p,q : polynomial degrees
%                   U,V : knot vectors in u,v-direction
%                    CP : Control Point coordinates and weights
%                   u,v : Parametric coordinate of the NURBS surface
%                     d : element displacement vector
%
%                Output :
% displacement_gradient : element strain vector strain vector epsilon = [exx eyy 2exy]'
%
% Function layout :
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

% Number of basis functions per element
ne = (p+1)*(q+1);  

% Number of degrees of freedom per element
ndof = 2*ne;

% 1. Compute the basis functions and its derivatives
[~,dR] = nurbs_basis_functions_and_first_derivatives_2D(i,p,u,U,j,q,v,V,CP);    

% 2. Compute the Jacobian of the transformation from the physical space to the
%    NURBS parameter domain:

% Initialize Jacobian
Jxxi = zeros(2,2);
% initialize counter
k = 0;
% Loop over all the non-zero contributions at the span under study
for c = 0:q
    for b = 0:p
        % Update counter
        k = k + 1;
        % Compute recursively the entries of the Jacobian
        Jxxi(1,1) = Jxxi(1,1) + CP(i-p+b,j-q+c,1)*dR(k,1);
        Jxxi(1,2) = Jxxi(1,2) + CP(i-p+b,j-q+c,2)*dR(k,1);
        Jxxi(2,1) = Jxxi(2,1) + CP(i-p+b,j-q+c,1)*dR(k,2);
        Jxxi(2,2) = Jxxi(2,2) + CP(i-p+b,j-q+c,2)*dR(k,2);
    end
end

% 3. Compute the derivatives of the NURBS basis functions on the physical
%    space from those on the parent domain

% Initialize matrix containing the derivatives for each basis function
dRdx = zeros(ne,2);

% Assign its entries recursively
for i=1:ne
    dRdx(i,:) = Jxxi\dR(i,:)';
end

% 4. Compute the B-operator matrices for both components of the displacement
% fields

% Initialize matrices
B_grad_u = zeros(2,ndof);
B_grad_v = zeros(2,ndof);

% Assign its entries recursively
for i=1:ne
    B_grad_u(1,2*i-1) = dRdx(i,1);
    B_grad_u(2,2*i-1) = dRdx(i,2);
    B_grad_v(1,2*i) = dRdx(i,1);
    B_grad_v(2,2*i) = dRdx(i,2);
end

% 5. Compute the displacement gradient matrix

displacement_gradient_u = B_grad_u*d;
displacement_gradient_v = B_grad_v*d;

% Initialize matrix
displacement_gradient = zeros(3,3);

% Assign its entries
displacement_gradient(1,1) = displacement_gradient_u(1);
displacement_gradient(1,2) = displacement_gradient_u(2);
displacement_gradient(2,1) = displacement_gradient_v(1);
displacement_gradient(2,2) = displacement_gradient_v(2);

end