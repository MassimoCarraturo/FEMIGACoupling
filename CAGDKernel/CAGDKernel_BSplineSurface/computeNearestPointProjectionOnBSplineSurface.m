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
function [xi,eta,Projected,flag,nIter] = computeNearestPointProjectionOnBSplineSurface(P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,newtonRaphson)
%% Function documentation
% 
% Returns the nearest point projection of a given point P onto a given
% NURBS-parametrized surface. The applied method for the solution of the
% non-linear system is the Newton-Rapson iterations.
%
% Source : Tianyang Wang (2010),'Investigation and Implementation of
%          Non-Matching Grids Data Transfer', Master Thesis, Technische
%          Universität München
%
%         Input :
%             P : The point to be projected on the B-Spline surface
%           p,q : The polynomial degrees of the NURBS surface
%        Xi,Eta : The knot vectors of the NURBS surface
%            CP : Set of Control Point coordinates and weights
%       isNURBS : Flag on whether the patch is a NURBS or a B-Spline
%      xi0,eta0 : Initial guess for the surface parameters
% newtonRaphson :   newtonRaphson.eps : Residual tolerance
%                 newtonRaphson.maxIt : Maximum number of iterations
%
%      Output :
%      xi,eta : The computed surface parameters for the closest point 
%               projection onto the NURBS surface
%   Projected : The closest point projection onto the NURBS surface
%        flag : 0 if the algorithm didnt converge or 1 if it converged
%       nIter : Number of iterations for convergence
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the Newton-Rapson iterations
%
%    1i. Update the iteration counter
%
%   1ii. Find the respective span
%
%  1iii. Compute the NURBS basis functions and their first and second derivatives
%
%   1iv. Compute the point on the NURBS surface
%
%    1v. Check for point coincidence
%
%   1vi. Compute the base vectors and their derivatives for the NURBS surface at (xi,eta)
%
%  1vii. Check the condition of zero cosine
%
% 1viii. Compute the updated values by applying the Newton's scheme
%
%   1ix. Modify (xi,eta) if they are out of their intervals of restriction
%
%    1x. Check condition for convergence
%
%% Function main bpdy

%% 0. Read input

% Initialize counter
counter = 0;

% Initialize the surface parameters
xi = xi0;
eta = eta0;

% Initialize flag to true
flag = 1;

% Compute the number of Control Points at each parametric direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

%% 1. Loop over all the Newton-Rapson iterations
while counter <= newtonRaphson.maxIt
    %% 1i. Update the iteration counter
    counter = counter + 1;
    
    %% 1ii. Find the respective span
    xiSpan = findKnotSpan(xi,Xi,nxi);
    etaSpan = findKnotSpan(eta,Eta,neta); 
    
    %% 1iii. Compute the NURBS basis functions and their first and second derivatives
    dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,2);
    
    %% 1iv. Compute the point on the NURBS surface
    Projected = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
    
    %% 1v. Check for point coincidence
    distance = (Projected - P);
    if norm(distance) <= newtonRaphson.eps
        break;
    end
    
    %% 1vi. Compute the base vectors and their derivatives for the NURBS surface at (xi,eta)
    [dG1,dG2] = computeBaseVectorsAndDerivativesForBSplineSurface(xiSpan,p,etaSpan,q,CP,1,dR);
    
    %% 1vii. Check the condition of zero cosine
    
    % Compute the cosine with respect to the u-parametric coordinate
    xiCos = abs(dG1(:,1)'*distance)/norm(dG1(:,1))/norm(Projected-P);
    
    % Compute the cosine with respect to the u-parametric coordinate
    etaCos = abs(dG2(:,1)'*distance)/norm(dG2(:,1))/norm(Projected-P);
    
    % Check the orthogonality condition
    if xiCos <= newtonRaphson.eps && etaCos <= newtonRaphson.eps
       break; 
    end
    
    %% 1viii. Compute the updated values by applying the Newton's scheme
    
    % Compute the system matrix
    J = [norm(dG1(:,1))^2 + distance'*dG1(:,2)    dG1(:,1)'*dG2(:,1) + distance'*dG1(:,3)
         dG2(:,1)'*dG1(:,1) + distance'*dG1(:,3)  norm(dG2(:,1))^2 + distance'*dG2(:,2)];
     
    % Compute the right-hand side residual
    r = - [distance'*dG1(:,1); distance'*dG2(:,1)];
    
    % Solve the linear system
    delta = J\r;
    
    % Update the parametric locations
    xi = xi + delta(1,1);
    eta = eta + delta(2,1);
    
    %% 1ix. Modify (xi,eta) if they are out of their intervals of restriction
    if xi>Xi(length(Xi)) 
        xi = Xi(length(Xi)) ; 
    end 
    if xi<Xi(1) 
        xi = Xi(1); 
    end
    if eta>Eta(length(Eta))
        eta = Eta(length(Eta)); 
    end
    if eta<Eta(1) 
        eta = Eta(1); 
    end
    
    %% 1x. Check additional condition for convergence
    condition = norm(delta(1,1)*dG1(:,1) + delta(2,1)*dG2(:,1));
    if condition <= newtonRaphson.eps
        break;
    end
end

%% 2. Assign the output values
nIter = counter - 1;
if nIter == newtonRaphson.maxIt
    flag = 0;
end

end