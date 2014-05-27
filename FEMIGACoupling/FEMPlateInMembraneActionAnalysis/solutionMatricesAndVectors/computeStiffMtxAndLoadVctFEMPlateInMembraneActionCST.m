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
function [K,F,minElSize] = computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST(u,uSaved,uDot,uDotSaved,...
    DOFNumbering,strMsh,analysis,F,bodyForces,strDynamics,parameters,int,outMsg)
%% Function documentation
%
% Returns the stiffness matrix and the load vector corresponding to the
% plate in membrane action analysis using the Constant Strain Triangle
% (CST) for the displacement field discretization.
%
%           Input :
%          uSaved : The discrete solution field of the previous 
%                   time step
%       uDotSaved : The time derivative of the discrete solution 
%                   field of the previous time step
%      uDDotSaved : The second order time derivative of the 
%                   discrete solution field of the previous time (dummy 
%                   variable for this function)
%          strMsh : Nodes and elements in the mesh
%        analysis : .type : The analysis type
%               F : The global load vector corresponding to surface
%                   tractions
%      bodyForces : Function handle to body force vector computation
%      parameters : Problem specific technical parameters
%             int : On the quadrature (numnerical integration)
%          outMsg : On outputting information
%
%          Output :
%               K : The master stiffness matrix of the system
%               F : The updated load vector accounting also for the body
%                   forces
%       minElSize : The minimum element area size in the mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Numnerical quadrature
%
% 2. Loop over all the elements in the mesh
% ->
%    2i. Get the element in the mesh
%
%   2ii. Get the nodes in the element
%
%  2iii. Create an Element Freedome Table (EFT)
%
%   2iv. Loop over all the quadrature points
%   ->
%        2iv.1. Transform the Gauss Point location from the parameter to the physical space
%
%        2iv.2. Compute the basis functions and their derivatives at the Gauss Point
%
%        2iv.3. Form the basis functions matrix at the Gauss Point
%
%        2iv.4. Form the B-Operator matrix for the plate in membrane action problem
%
%        2iv.5. Compute the element stiffness matrix at the Gauss Point and assemble to master stiffness matrix via the EFT
%
%        2iv.6. Compute the element load vector due to body forces and assemble to master load vector via the EFT
%   <-
% <-
%
% 3. Appendix
%
%% Functions main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________\n');
    fprintf('###################################################################\n');
    fprintf('Computation the master stiffness matrix for a FEM discretized plate\n');
    fprintf('in membrane action problem has been initiated\n\n');
    fprintf('Analysis type : plane ');
    if strcmp(analysis.type,'PLANE_STRESS')
        fprintf('stress\n');
    elseif strcmp(analysis.type,'PLANE_STRAIN')
        fprintf('strain\n');
    end
    fprintf('___________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of nodes in the mesh
nNodes = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
nDOFs = 2*nNodes;

% Number of nodes at the element level
nNodesEl = 3;

% Number of DOFs at the element level
nDOFsEl = 2*nNodesEl;

% Minimum element edge size
minElSize = 'undefined';

% Compute the material matrix for the given problem
if strcmp(analysis.type,'PLANE_STRESS')
    preFactor = parameters.E/(1-parameters.nue^2);
    C = preFactor*[1             parameters.nue 0
                   parameters.nue 1              0
                   0              0             (1-parameters.nue)/2];
elseif strcmp(analysis.physics,'PLANE_STRAIN')
    preFactor = parameters.E*(1-parameters.nue)/(1+parameters.nue)/(1-2*parameters.nue);
    C = preFactor*[1                                 parameters.nue/(1-parameters.nue) 0
                   parameters.nue/(1-parameters.nue) 1                                 0
                   0                                 0                                 (1-2*parameters.nue)/2/(1-parameters.nue)];
end
    
% Initialize output array
K = zeros(nDOFs,nDOFs);

%% 1. Numnerical quadrature
if strcmp(int.type,'default')
    nGP = 1;
elseif strcmp(int.type,'manual')
    nGP = int.nGP;
end
[GP,GW] = getGaussRuleOnCanonicalTriangle(nGP);

%% 2. Loop over all the elements in the mesh
for counterEl = 1:length(strMsh.elements(:,1))
    %% 2i. Get the element in the mesh
    element = strMsh.elements(counterEl,:);
    
    %% 2ii. Get the nodes in the element
    node1 = strMsh.nodes(element(1,1),:);
    node2 = strMsh.nodes(element(1,2),:);
    node3 = strMsh.nodes(element(1,3),:);
    
    %% 2iii. Create an Element Freedome Table (EFT)
    EFT = zeros(nDOFsEl,1);
    for counterEFT = 1:nNodesEl
        EFT(2*counterEFT-1) = 2*element(1,counterEFT)-1;
        EFT(2*counterEFT) = 2*element(1,counterEFT);
    end
    
    %% 2iv. Loop over all the quadrature points
    for counterGP = 1:nGP
        %% 2iv.1. Transform the Gauss Point location from the parameter to the physical space
        xGP = GP(counterGP,1)*node1(1,:) + GP(counterGP,2)*node2(1,:) + (1-GP(counterGP,1)-GP(counterGP,2))*node3(1,:);
        
        %% 2iv.2. Compute the basis functions and their derivatives at the Gauss Point
        [dN,detJxxi] = computeCST2DBasisFunctionsAndFirstDerivatives(node1,node2,node3,xGP(1,1),xGP(1,2));
        
        %% 2iv.3. Form the basis functions matrix at the Gauss Point
        N = [dN(1,1) 0       dN(2,1) 0       dN(3,1) 0
             0       dN(1,1) 0       dN(2,1) 0       dN(3,1)];
        
        %% 2iv.4. Form the B-Operator matrix for the plate in membrane action problem
        B = [dN(1,2) 0       dN(2,2) 0       dN(3,2) 0
             0       dN(1,3) 0       dN(2,3) 0       dN(3,3)
             dN(1,3) dN(1,2) dN(2,3) dN(2,2) dN(3,3) dN(3,2)];
         
        %% 2iv.5. Compute the element stiffness matrix at the Gauss Point and assemble to master stiffness matrix via the EFT
        K(EFT,EFT) = K(EFT,EFT) + (B'*C*B)*detJxxi*GW(counterGP);
        
        %% 2iv.6. Compute the element load vector due to body forces and assemble to master load vector via the EFT
        bF = bodyForces(xGP(1,1),xGP(1,2),xGP(1,3));
        F(EFT) = F(EFT) + N'*bF(1:2,1)*detJxxi*GW(counterGP);
    end
end

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nComutation of the master stiffness matrix took %.2d seconds \n\n',computationalTime);
    fprintf('_________________Stiffness Matrix Computation Ended_______________\n');
    fprintf('##################################################################\n\n\n');
end

end