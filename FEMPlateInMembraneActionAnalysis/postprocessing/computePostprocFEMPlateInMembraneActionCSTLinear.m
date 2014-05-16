%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function [epsilon,sigma] = computePostprocFEMPlateInMembraneActionCSTLinear(strMsh,analysis,parameters,dHat)
%% Function documentation
%
% Returns the strain field [epsilonXX epsilonYY epsilonXY] and stress field
% [sigmaXX sigmaYY sigmaXY] corresponding to a linear classical finite 
% element plate in membrane action analysis using the constant strain
% triangle.
%
%       Input :
%      strMsh : Nodes and elements in the mesh
%    analysis : .type : The analysis type
%  parameters : Problem specific technical parameters
%        dHat : The nodal displacement field
%
%      Output :
%     epsilon : The strain field [epsilonXX epsilonYY epsilonXY]
%       sigma : The stress field [sigmaXX sigmaYY sigmaXY]
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the elements in the mesh
% ->
%    1i. Get the element in the mesh
%
%   1ii. Get the nodes in the element
%
%  1iii. Create an Element Freedome Table (EFT)
%
%   1iv. Compute the basis functions and their derivatives at the Gauss Point
%
%    1v. Form the B-Operator matrix for the plate in membrane action problem
%
%   1vi. Compute the strain vector in a Voigt notation at the current element
%
%  1vii. Compute the stress vector in a Voigt notation at the current element
%
% 1viii. Correct the shear component of the strain epsilon = [epsilonXX epsilonYY 2*epsilonXY]
% <-
%
%% Function main body

%% 0. Read input

% Number of nodes at the element level
nNodesEl = 3;

% Number of DOFs at the element level
nDOFsEl = 2*nNodesEl;

% Number of nodes in the mesh
noElements = length(strMsh.elements);

% Compute the material matrix for the given problem (The shear entry is 
% multiplied by two so that it returns the true strain field and not the 
% one needed for the computation of the internal virtual work)
if strcmp(analysis.type,'PLANE_STRESS')
    preFactor = parameters.E/(1-parameters.nue^2);
    C = preFactor*[1             parameters.nue  0
                   parameters.nue 1              0
                   0              0             (1-parameters.nue)/2];
elseif strcmp(analysis.physics,'PLANE_STRAIN')
    preFactor = materialProperties.E*(1-parameters.nue)/(1+parameters.nue)/(1-2*parameters.nue);
    C = preFactor*[1                                 parameters.nue/(1-parameters.nue) 0
                   parameters.nue/(1-parameters.nue) 1                                 0
                   0                                 0                                 (1-2*parameters.nue)/2/(1-parameters.nue)];
end

% Initialize output arrays
epsilon = zeros(3,noElements);
sigma = zeros(3,noElements);

%% 1. Loop over all the elements in the mesh
for counterEl = 1:length(strMsh.elements(:,1))
    %% 1i. Get the element in the mesh
    element = strMsh.elements(counterEl,:);
    
    %% 1ii. Get the nodes in the element
    node1 = strMsh.nodes(element(1,1),:);
    node2 = strMsh.nodes(element(1,2),:);
    node3 = strMsh.nodes(element(1,3),:);
    
    %% 1iii. Create an Element Freedome Table (EFT)
    EFT = zeros(nDOFsEl,1);
    for counterEFT = 1:nNodesEl
        EFT(2*counterEFT-1) = 2*element(1,counterEFT)-1;
        EFT(2*counterEFT) = 2*element(1,counterEFT);
    end
    
    %% 1iv. Compute the basis functions and their derivatives at the Gauss Point
	[dN,~] = computeCST2DBasisFunctionsFirstDerivatives(node1,node2,node3);
    
    %% 1v. Form the B-Operator matrix for the plate in membrane action problem
    B = [dN(1,1)    0        dN(2,1)   0        dN(3,1)    0
         0          dN(1,2)  0         dN(2,2)  0          dN(3,2)
         dN(1,2)    dN(1,1)  dN(2,2)   dN(2,1)  dN(3,2)    dN(3,1)];
     
    %% 1vi. Compute the strain vector in a Voigt notation at the current element
    epsilon(:,counterEl) = B*dHat(EFT);
    
    %% 1vii. Compute the stress vector in a Voigt notation at the current element
    sigma(:,counterEl) = C*epsilon(:,counterEl);
    
    %% 1viii. Correct the shear component of the strain epsilon = [epsilonXX epsilonYY 2*epsilonXY]
    epsilon(3,counterEl) = epsilon(3,counterEl)/2;
end

end