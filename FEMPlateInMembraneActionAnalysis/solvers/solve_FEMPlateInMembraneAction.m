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
function [dHat,FComplete,minElSize] = ...
    solve_FEMPlateInMembraneAction(strMsh,homDOFs,inhomDOFs,valuesInhomDOFs,...
    NBC,F,bodyForces,analysis,parameters,nLinearAnalysis,strDynamics,int,...
    caseName,pathToOutput,outMsg)
%% Function documentation
%
% Returns the displacement field, the complete force vector and the minimum
% element area size for a plate in membrane action problem solved with the
% classical Finite Element Method.
%
%           Input :
%          strMsh : Nodes and elements in the mesh
%         homDOFs : The global numbering of the nodes where homogeneous
%                   Dirichlet boundary conditions are applied 
%       inhomDOFs : The global numbering of the nodes where inhomogeneous
%                   Dirichlet boundary conditions are applied
% valuesInhomDOFs : Prescribed values on the nodes where inhomogeneous
%                   Dirichlet boundary conditions are applied
%             NBC :    .nodes : The nodes where Neumann boundary 
%                                conditions are applied
%                    .loadType : The type of the load for each Neumann node
%                   .fctHandle : The function handle for each Neumann node
%                                for the computation of the load vector 
%                                (these functions are unde the folder load)
%               F : Global load vector corresponding to surface tractions
%      bodyForces : Function handle to body force vector computation
%        analysis : .type : The analysis type
%      parameters : Problem specific technical parameters
% nLinearAnalysis :     .scheme : The employed nonlinear scheme
%                    .tolerance : The residual tolerance
%                      .maxIter : The maximum number of the nonlinear 
%                                 iterations
%     strDynamics : .timeDependence : Steady-state or transient analysis
%                           .scheme : The time integration scheme
%                               .T0 : The start time of the simulation
%                             .TEnd : The end time of the simulation
%                          .nTSteps : The number of the time steps
%             int : On the numerical integration (quadrature)
%                      .type : 'default', 'manual'
%                       .nGP : Number of Gauss Points
%        caseName : The name of the case in the inputGiD case folder
%    pathToOutput : Define the path to where to write out the results
%          outMsg : On outputting information
%
%          Output :
%            dHat : The nodal displacement field
%       FComplete : The complete force vector
%       minElSize : The minimum element area size in the mesh
%
% Function layout :
%
% 
%
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________\n');
    fprintf('#############################################################\n');
    fprintf('Compute the displacement field for a plate in membrane action\n');
    fprintf('problem has been initiated\n');
    fprintf('_____________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of nodes in the mesh
noNodes = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
nDOFs = 2*noNodes;

% GLobal DOF numbering
DOFNumbering = 1:nDOFs;

% Assign dummy variables
uSaved = 'undefined';
uDot = 'undefined';
uDotSaved = 'undefined';
uDDotSaved = 'undefined';

% Initialize output array
dHat = zeros(nDOFs,1);

% Title for the output file
title = 'geometrically linear steady-state plane stress analysis';

% Connectivity arrays for the DOFs into the resulting vectors
xDisp = 1:2:nDOFs-1;
yDisp = 2:2:nDOFs;

% Make directory to write out the results of the analysis
isExistent = exist(strcat('../../outputVTK/FEMPlateInMembraneActionAnalysis/',caseName),'dir');
if ~isExistent
    mkdir(strcat('../../outputVTK/FEMPlateInMembraneActionAnalysis/',caseName));
end

%% 1. Find the prescribed and the free DOFs of the system

% Prescribed DOFs (DOFs on which either homogeneous or inhomogeneous 
% Dirichlet boundary conditions are prescribed)
prescribedDoFs = mergesorted(homDOFs,inhomDOFs);
prescribedDoFs = unique(prescribedDoFs);

% Free DOFs of the system (actual DOFs over which the solution is computed)
freeDOFs = DOFNumbering;
freeDOFs(ismember(freeDOFs,prescribedDoFs)) = [];

%% 2. Solve the linear equation system
[dHat,~,~,FComplete,minElSize] = solve_FEMLinearSystem(uSaved, ...
    uDotSaved,uDDotSaved,strMsh,analysis,F,bodyForces,parameters,dHat,uDot, ...
    @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST,DOFNumbering,...
    freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,strDynamics,nLinearAnalysis,int,'');

%% 3. Postprocessing
[epsilon,sigma] = computePostprocFEMPlateInMembraneActionCSTLinear(strMsh,analysis,parameters,dHat);

%% 4. Write out the results into a file
fprintf('>> Writting out the results to "%s"\n',strcat(pathToOutput,caseName,'/'));
nodalDisplacement = [dHat(xDisp)'; dHat(yDisp)'; zeros(1,noNodes)];
writeOutputFEMPlateInMembraneActionToVTK(strMsh,nodalDisplacement,epsilon,sigma,caseName,pathToOutput,title);

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nLinear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('_____________________Linear Analysis Ended___________________\n');
    fprintf('#############################################################\n\n\n');
end

end