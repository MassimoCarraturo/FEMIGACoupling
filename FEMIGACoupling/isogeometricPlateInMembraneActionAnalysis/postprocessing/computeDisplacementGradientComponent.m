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
function displacement_gradient= get_displacement_gradient_component(i,p,u,U,j,q,v,V,CP,d,error)
%% Function documentation
%
% Returns the displacement gradient vectors = [ui,x ui,y]' at
% the parametric location (u,v), i=x,y
%
%                 Input : 
%                   p,q : polynomial degrees
%                   U,V : knot vectors in u,v-direction
%                    CP : Control Point coordinates and weights
%                   u,v : Parametric coordinate of the NURBS surface
%                     d : element displacement vector
%                 error : structure containing all information on the error
%                         computation
%
%                Output :
% displacement_gradient : element displacement gradient tensor
%
% Function layout :
%
% 1. Compute the basis functions and its derivatives
%
% 2. Compute the Jacobian of the transformation from the physical space to the
%    NURBS parameter domain:
%
% 3. Compute the B-operator matrix
%
% 4. Compute the displacement gradient component
%
%% Function main body

% Number of basis functions per element
ne = (p+1)*(q+1);  

% Number of degrees of freedom per element
ndof = 2*ne;

% 1. Compute the basis functions and its derivatives
[~,dR] = nurbs_basis_functions_and_first_derivatives_2D(i,p,u,U,j,q,v,V,CP);    

% 2. Compute the Jacobian of the transformation from the physical space to the
% NURBS parameter domain:

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

% 4. Compute the B-operator matrix
% initialize matrix
B_grad = zeros(2,ndof);

% compute its entries recursively
for i=1:ne
    B_grad(1,2*i+error.component-2) = dRdx(i);
    B_grad(2,2*i+error.component-2) = dRdx(i);
end

% 5. Compute the displacement gradient component
displacement_gradient = B_grad*d;

end