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
function [Ke,Fe] = computeElStiffMtxAndLoadVctPlateInMembraneActionLinear(noCPsLoc,R,dR,Jxxi,D,FBody)
%% Function documentation
%
% Evaluates the stiffness matrix corresponding to the isogeometric plate in 
% membrane action of an element at the quadrature point
% 
%      Input :
%   noCPsLoc : Local number of Control Points
%          R : The NURBS basis functions
%         dR : The vector of the derivatives of the NURBS basis functions
%              at the quadrature point
%       Jxxi : jacobian of the transformation from the physical space the 
%              NURBS parametric space evaluated at the quadrature point
%          D : material matrix
%      FBody : The body force vector in the parametric location where the
%              basis functions where computed
%
%     Output :
%         Ke : The element stiffness matrix
%         Fe : The element load vector due to the body forces
%
% Function layout :
%
% 1. Compute the body force vector at the parametric location where the basis functions are computed
%
% 2. Compute the derivatives of the basis functions w.r.t. the physical domain at the quadrature point
%
% 3. Compute B-operator matrix at the parametric location where the basis functions are computed
%
% 4. Compute the element stiffness matrix at the parametric location where the basis functions are computed
%
%% Function main body

%% 1. Compute the body force vector at the parametric location where the basis functions are computed
if norm(FBody)~=0
   RMatrix = zeros(2,2*noCPsLoc);
   for i = 1:noCPsLoc
        RMatrix(1,2*i-1) = R(i,1);
        RMatrix(2,2*i) = R(i,1);
   end
   Fe = RMatrix*FBody;
else
   Fe = zeros(2*noCPsLoc,1); 
end

%% 2. Compute the derivatives of the basis functions w.r.t. the physical domain at the quadrature point

% Initialize matrix containing the derivatives for each basis function
dRdx = zeros(noCPsLoc,2);

% Compute the derivatives on the physical space given those at the parent
% domain
for i=1:noCPsLoc
    dRdx(i,:) = Jxxi\dR(i,:)';
end

%% 3. Compute B-operator matrix at the parametric location where the basis functions are computed

% Initialize B-operator matrix
B = zeros(3,2*noCPsLoc);

% Loop over all the entries
for i=1:noCPsLoc
   B(1,2*i-1) = dRdx(i,1);
   B(2,2*i) = dRdx(i,2);
   B(3,2*i-1) = dRdx(i,2);
   B(3,2*i) = dRdx(i,1);
end

%% 4. Compute the element stiffness matrix at the parametric location where the basis functions are computed
Ke = B'*D*B;

end