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
function plot_postprocResultantsIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,parameters,xiGrid,etaGrid,dHat,graph)
%% Function documentation
%
% Plots the reference configuration of an isogeometric plate in membrane 
% action together with the selected postprocessing resultant.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
%     parameters : Technical and geometrical parameters of the plate
%           dHat : The displacement field of the control points
% xiGrid,etaGrid : The grid points used for the plotting of the NURBS
%                  geometry
%           dHat : The displacement field of the control points
%          graph : Information on the graphics
%
%         Output :
%                  graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the resultant array to be used for the visualization
%
% 2. Visualize the resultant and the knots over the domain
%
%% Function main body

%% 0. Read input

% Number of knots in u,v-direction
mxi = length(Xi);
meta = length(Eta);

% Assign a tolerance value
tol = 10e-10;

% Compute incremental steps for the resultant computation over the domain
% incremental step for eta:
deta = (Eta(meta)-Eta(1))/(xiGrid-1);

% incremental step for xi:
dxi = (Xi(mxi)-Xi(1))/(etaGrid-1);

% Initialize array of the resultant to be visualized
resultantComponent = zeros(xiGrid,etaGrid);

% Cartesian image of the parameter space
P = zeros(xiGrid,etaGrid,3);

% Number of Control Points in xi-,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Compute the material matrix for the plane stress problem
D = parameters.E/(1-parameters.nue^2)*[1 parameters.nue 0; parameters.nue 1 0; 0 0 (1-parameters.nue)/2];

% Make a DOF numbering for the given patch
dHatElem = zeros(mxi-p-1,meta-q-1,2*(p+1)*(q+1));
for etaSpan = (q+1):(meta-q-1)
    for xiSpan = (p+1):(mxi-p-1)
        xiCounter = 1; 
        for c = etaSpan-q-1:etaSpan-1 
            for b = xiSpan-p:xiSpan
                dHatElem(xiSpan,etaSpan,xiCounter)   = dHat(2*(c*nxi+b)-1);
                dHatElem(xiSpan,etaSpan,xiCounter+1) = dHat(2*(c*nxi+b));
                xiCounter = xiCounter + 2;
            end
        end
    end
end

% Initialize flags
isUndeformed = 0;
isDeformed = 0;

%% 1. Compute the resultant array to be used for the visualization

% counting index in eta-direction
etaCounter = 1;  

% Initialize coordinate in eta-direction
eta = Eta(1);

% Loop over all the parametric coordinates in eta-direction
while eta <= Eta(meta)+tol
    % Find the span in the eta-direction
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    % Initialize coordinate in xi-direction
    xi = Xi(1);
    
    % Initialize counter in xi-direction
    xiCounter = 1;
    
    % Loop over all the parametric coordinates in xi-direction
    while xi <= Xi(mxi)+tol
        % Find the span in xi-direction
        xiSpan = findKnotSpan(xi,Xi,nxi);
        
        % Compute the IGA basis functions and possibly their derivatives
        if strcmp(graph.resultant,'displacement')
            nDrv = 0;
        else
            nDrv = 1;
        end
        dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
        
        % Compute the Cartesian image of the paratric point
        P(xiCounter,etaCounter,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
        
        % Get the actual Control Point displacement vector for the current
        % knot span
        dHatActual(:,1) = dHatElem(xiSpan,etaSpan,:);
        
        % Compute the resultant component
        if strcmp(graph.resultant,'displacement')
            % Compute the displacement field at the given parametric
            % location
            d = computePostprocDisplacementIGAPlateInMembraneAction(p,q,dR(:,1),dHatActual);
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'x')
                resultantComponent(xiCounter,etaCounter) = d(1,1);
            elseif strcmp(graph.component,'y')
                resultantComponent(xiCounter,etaCounter) = d(2,1);
            elseif strcmp(graph.component,'2norm')
                resultantComponent(xiCounter,etaCounter) = sqrt(d(1,1)^2 + d(2,1)^2);
            end
        elseif strcmp(graph.resultant,'strain')
            % Compute the strain tensor in a Voigt notation
            epsilonVoigt = computePostprocVoigtStrainIGAPlateInMembraneAction(xiSpan,p,etaSpan,q,CP,dR(:,2:3),dHatActual);
             
            % Compute the shear strain component, recall 
            % epsilon = [epsilonXX epsilonYY 2*epsilonXY]' 
            epsilonVoigt(3,1) = epsilonVoigt(3,1)/2;
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'x')
                resultantComponent(xiCounter,etaCounter) = epsilonVoigt(1,1);
            elseif strcmp(graph.component,'y')
                resultantComponent(xiCounter,etaCounter) = epsilonVoigt(2,1);
            elseif strcmp(graph.component,'xy')
                resultantComponent(xiCounter,etaCounter) = epsilonVoigt(3,1);
            elseif strcmp(graph.component,'1Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(epsilonVoigt(1,1)+epsilonVoigt(2,1) + sqrt((epsilonVoigt(1,1)-epsilonVoigt(2,1))^2 + epsilonVoigt(3,1)^2));
            elseif strcmp(graph.component,'2Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(epsilonVoigt(1,1)+epsilonVoigt(2,1) - sqrt((epsilonVoigt(1,1)-epsilonVoigt(2,1))^2 + epsilonVoigt(3,1)^2));
            end
        elseif strcmp(graph.resultant,'stress')
            % Compute the strain tensor in a Voigt notation
            epsilonVoigt = computePostprocVoigtStrainIGAPlateInMembraneAction(xiSpan,p,etaSpan,q,CP,dR(:,2:3),dHatActual);
            
            % Compute the stress tensor in a Voigt notation
            sigmaVoigt = D*epsilonVoigt;
            
            % Decide upon the the visualization resultant
            if strcmp(graph.component,'x')
                resultantComponent(xiCounter,etaCounter) = sigmaVoigt(1,1);
            elseif strcmp(graph.component,'y')
                resultantComponent(xiCounter,etaCounter) = sigmaVoigt(2,1);
            elseif strcmp(graph.component,'xy')
                resultantComponent(xiCounter,etaCounter) = sigmaVoigt(3,1);
            elseif strcmp(graph.component,'1Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(sigmaVoigt(1,1)+sigmaVoigt(2,1) + sqrt((sigmaVoigt(1,1)-sigmaVoigt(2,1))^2 + 4*sigmaVoigt(3,1)^2));
            elseif strcmp(graph.component,'2Principal')
                resultantComponent(xiCounter,etaCounter) = 0.5*(sigmaVoigt(1,1)+sigmaVoigt(2,1) - sqrt((sigmaVoigt(1,1)-sigmaVoigt(2,1))^2 + 4*sigmaVoigt(3,1)^2));
            end
        end

        % Update the counter in xi-direction
        xiCounter = xiCounter + 1;
        
        % Update the parametric coordinate in xi-direction
        xi = xi + dxi;
    end
    % Update the counter in eta-direction
    etaCounter = etaCounter + 1;
    
    % Update the parametric coordinate in eta-direction
    eta = eta + deta;
end

%% 2. Visualize the resultant and the knots over the domain

% Plot the resultant over the reference configuration
surf(P(:,:,1),P(:,:,2),P(:,:,3),resultantComponent(:,:));
hold on;

% Plot the element boundaries on the undeformed geometry
plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isUndeformed,xiGrid,etaGrid);
hold off;

end

