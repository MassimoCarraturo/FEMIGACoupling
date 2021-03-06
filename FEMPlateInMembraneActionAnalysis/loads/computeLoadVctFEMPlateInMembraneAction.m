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
function F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,int,outMsg)
%% Function documentation
%
% Returns the global load vector corresponding to classical Finite Element
% analysis for the case of a plate in membrane action problem.
%
%           Input :
%          strMsh : Nodes and elements in the mesh
%        analysis : .type : The analysis type
%             NBC :    .nodes : The nodes where Neumann boundary 
%                                conditions are applied
%                    .loadType : The type of the load for each Neumann node
%                   .fctHandle : The function handle for each Neumann node
%                                for the computation of the load vector 
%                                (these functions are unde the folder load)
%               t : The time instance
%             int : On the numerical quadrature scheme
%                       int.domain : 'straightLine' or 'triangular'
%                         int.type : 'default' or 'manual'
%                          int.nGP : The number of Gauss points for the
%                                    integration
%          outMsg : On outputting information
%
%          Output :
%               F : The global load vector
%
% Function layout :
%
% 0. Read input
%
% 1. Get the Gauss integration rules
%
% 2. Loop over all elements on the Neumann boundary
% ->
%    2i. Get the nodes which are on the Neumann boundary
%
%   2ii. Get the nodes of the element on the Neumann boundary
%
%  2iii. Get the function handle for this type of loading
%
%   2iv. Assemble the Element Freedome Table (EFT)
%
%    2v. Loop over all Gauss Points
%    ->
%        2v.1. Transform the Gauss Point from the parameter to the physical space
%
%        2v.2. Compute the basis functions at the Gauss Point
%
%        2v.3. Sort the basis functions into the basis functions array
%
%        2v.4. Compute the applied load onto the Gauss Point
%
%        2v.5. Compute the determinant of the transformation from the physical space to the parameter space
%
%        2v.6. Compute the element load vector on the Gauss point
%
%        2v.7. Assemble the contribution to the global load vector via the EFT
%    <-
% <-
% 3. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('______________________________________________________________\n');
    fprintf('##############################################################\n');
    fprintf('Computation the global load vector for a FEM discretized plate\n');
    fprintf('in membrane action problem has been initiated\n');
    fprintf('______________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of nodes
nNodes = length(strMsh.nodes(:,1));

% Number of DOFs
if strcmp(analysis.type,'PLANE_STRESS')
    nDOFs = 2*nNodes;
end

% Number of element nodes
nNodesEl = length(strMsh.nodes(1,:));

% Number of element DOFs
if strcmp(analysis.type,'PLANE_STRESS')
    nDOFsEl = 2*nNodesEl;
end

% Initialize output array
F = zeros(nDOFs,1);

%% 1. Get the Gauss integration rules

% For the boundary integration
if strcmp(int.type,'default')
    nGP = 1;
else
    nGP = int.nGP;
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(nGP);

%% 2. Loop over all elements on the Neumann boundary
for counterEl = 1:length(NBC.lines(:,1))
    %% 2i. Get the nodes which are on the Neumann boundary
    nodeIDs = NBC.lines(counterEl,1:2);
    x1 = strMsh.nodes(nodeIDs(1),:);
    x2 = strMsh.nodes(nodeIDs(2),:);
    
    %% 2ii. Get the nodes of the element on the Neumann boundary
    elementID = NBC.lines(counterEl,3);
    element = strMsh.elements(elementID,:);
    node1 = strMsh.nodes(element(1),:);
    node2 = strMsh.nodes(element(2),:);
    node3 = strMsh.nodes(element(3),:);
    
    %% 2iii. Get the function handle for this type of loading
    loadFctHandle = str2func(NBC.fctHandle(counterEl,:));
    
    %% 2iv. Assemble the Element Freedome Table (EFT)
    EFT = zeros(nDOFsEl,1);
    for counterEFT = 1: nNodesEl
        EFT(2*counterEFT-1) = 2*element(counterEFT)-1;
        EFT(2*counterEFT) = 2*element(counterEFT);
    end
    
    %% 2v. Loop over all Gauss Points
    for counterGP = 1:nGP
        %% 2v.1. Transform the Gauss Point from the parameter to the physical space
        xGP = (1 - GP(counterGP))*x1/2 + (1 + GP(counterGP))*x2/2;
        
        %% 2v.2. Compute the basis functions at the Gauss Point
        basisFctOnGP = computeCST2DBasisFunctions(node1,node2,node3,xGP(1,1),xGP(1,2));
        
        %% 2v.3. Sort the basis functions into the basis functions array
        N = [basisFctOnGP(1,1)  0                  basisFctOnGP(1,2) 0                 basisFctOnGP(1,3)    0
             0                  basisFctOnGP(1,1)  0                 basisFctOnGP(1,2) 0                    basisFctOnGP(1,3)];
        
        %% 2v.4. Compute the applied load onto the Gauss Point
        tractionOnGP = loadFctHandle(xGP(1,1),xGP(1,2),xGP(1,3),t);
        tractionOnGP2D = tractionOnGP(1:2,1);
        
        %% 2v.5. Compute the determinant of the transformation from the physical space to the parameter space
        
        % The linear mapping from the parameter to the physical space is
        % x(xi) = .5*(x2 - x1) + .5*(x2 + x1) => detJxxi = .5*(x2 - x1)
        detJxxi = norm(x2 - x1)/2;
        
        %% 2v.6. Compute the element load vector on the Gauss point
        FElOnGP = N'*tractionOnGP2D*detJxxi*GW(counterGP);
        
        %% 2v.7. Assemble the contribution to the global load vector via the EFT
        F(EFT) = F(EFT) + FElOnGP; 
    end
end

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nComutation of the load vector took %.2d seconds \n\n',computationalTime);
    fprintf('_________________Load Vector Computation Ended________________\n');
    fprintf('##############################################################\n\n\n');
end

end