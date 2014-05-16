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
function [u,flag] = computePointCoordinatesOnCanonicalBilinearQuadrilateral(x,xCorner,maxIt,eps)
%% Function documentation
%
% Returns the coordinates of a point with coordinates x within the 2D 
% quadrilateral x1-x2-x3-x4 in an arbitrary Cartesian space to the
% canonical bilinear quadrilateral. The applied method for solving the
% non-linear system is Newton-Rapson.
%
%       Input :
%           x : The coordinates of the point into the Cartesian system
%     xCorner : The coordinates of the corner nodes of the quadrilateral
%       maxIt : Maximum iterations for the Newton Rapson algorithm
%         eps : Tolerance for the Newton-Rapson convergence criterion based
%               on the 2-norm of the residual
%
%      Output :
%           u : The u,v coordinates of the point in the canonical space
%        flag : Boolean on the convergence of the Newton iterations
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the Newton iterations
%
%    1i. Compute the bilinear basis functions and its derivatives at u
%
%   1ii. Compute the right-hand side residual
%
%  1iii. Compute the Jacobian
%
%   1iv. Solve the linear system
%
%    1v. Check condition for convergence
%
%   1vi. Update the iteration counter
%
% 2. Check if convergence has been achieved
%
%% Function main body

%% 0. Read input

% Coordinates of point x into the canonical space. Initialize them into the
% center of the canonical bilinear quadrilateral
u = zeros(2,1);

% Initialize counter
counter = -1;

% Initialize convergence flag to true
flag = 1;

%% 1. Loop over all the Newton iterations
while counter<=maxIt
    %% 1i. Compute the bilinear basis functions and its derivatives at u
    [N,dN] = computeBilinearBasisFunctionsAndFirstDerivatives(u(1,1),u(2,1));
    
    %% 1ii. Compute the right-hand side residual
    
    % Compute the physical coordinates of the point with canonical
    % coordinates u
    xPredicted = zeros(2,1);
    DxPredictedDxi = zeros(2,1);
    DxPredictedDeta = zeros(2,1); 
    
    % Loop over all the basis functions and add contributions
    for i=1:length(N(1,:))
        % The coordinates of the point
        xPredicted = xPredicted + N(1,i)*xCorner(:,i);
        
        % The derivatives of the position vector w.r.t. tp xi
        DxPredictedDxi = DxPredictedDxi + dN(1,i)*xCorner(:,i);
        
        % The derivatives of the position vector w.r.t. tp eta
        DxPredictedDeta = DxPredictedDeta + dN(2,i)*xCorner(:,i);
    end
    
    % The residual
    residual = x - xPredicted;
    
    %% 1iii. Compute the Jacobian
    J = [DxPredictedDxi(1,1) DxPredictedDeta(1,1)
         DxPredictedDxi(2,1) DxPredictedDeta(2,1)];
    
    %% 1iv. Solve the linear system
    delta = J\residual;
    
    % Update the parametric locations
    u = u + delta;
    
    %% 1v. Check condition for convergence
    condition = norm(residual);
    
    if condition<=eps
        break;
    end
    
    %% 1vi. Update the iteration counter
    counter = counter + 1;
end

%% 2. Check if convergence has been achieved
if counter==maxIt
    flag = 0;
end

end

