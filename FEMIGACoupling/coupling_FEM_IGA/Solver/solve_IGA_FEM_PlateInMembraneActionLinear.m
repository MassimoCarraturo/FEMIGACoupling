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
function [dHat,FComplete,minElSize] = solve_IGA_FEM_PlateInMembraneActionLinear...
    (p,Xi,q,Eta,CP,isNURBS,parameters,Fl,bodyForces,homDOFs,inhomDOFs,valuesInhomDOFs,int,outMsg,...
    strMsh,homDOFsFEM,inhomDOFsFEM,valuesInhomDOFsFEM,...
    NBC,F,bodyForcesFEM,analysis,parametersFEM,nLinearAnalysisFEM,strDynamics,intFEM,...
    caseName,pathToOutput,outMsgFEM);

%% Function documentation
%
% Returns the displacement field  and the complete force vector 
% corresponding to 2D isogeometric plate in membrane action
%
%           Input :
%             p,q : polynomial degrees
%          Xi,Eta : knot vectors
%              CP : vector containing control point coordinates and weights
%         isNURBS : Flag on whether the geometrical basis is NURBS or 
%                   B-Spline
%      parameters : Technical and geometrical parameters
%              Fl : force vector
%      bodyForces : Handle to function for the computation of the body 
%                   forces
%         homDOFs : Array containing information on homogeneous Dirichlet
%                   boundary conditions
%   	inhomDOFs : Array containing information on the inhomogeneous 
%                   Dirichlet boundary conditions
% valuesInhomDOFs : Values on the inhomogeneous Dirichlet boundary 
%                   conditions
%             int : On the numerical integration
%          outMsg : Whether or not to output message on refinement progress
%                   'outputEnabled' : enables output information
%
%          Output :
%            dHatFEM : the displacement field
%       FComplete : the complete load vector
%       minElSize : The minimum element area size in the isogeometric mesh
%          strMsh : Nodes and elements in the mesh
%         homDOFsFEM : The global numbering of the nodes where homogeneous
%                   Dirichlet boundary conditions are applied 
%       inhomDOFsFEM : The global numbering of the nodes where inhomogeneous
%                   Dirichlet boundary conditions are applied
% valuesInhomDOFsFEM : Prescribed values on the nodes where inhomogeneous
%                   Dirichlet boundary conditions are applied
%             NBC :    .nodes : The nodes where Neumann boundary 
%                                conditions are applied
%                    .loadType : The type of the load for each Neumann node
%                   .fctHandle : The function handle for each Neumann node
%                                for the computation of the load vector 
%                                (these functions are unde the folder load)
%             IBC :     .nodes : The nodes where boundary 
%                                conditions are applied
%                    .boundaryType : The type of the boundary for each node
%                   
%               F : Global load vector corresponding to surface tractions
%      bodyForcesFEM : Function handle to body force vector computation
%        analysis : .type : The analysis type
%      parametersFEM : Problem specific technical parameters
% nLinearAnalysisFEM :     .scheme : The employed nonlinear scheme
%                    .tolerance : The residual tolerance
%                      .maxIter : The maximum number of the nonlinear 
%                                 iterations
%     strDynamics : .timeDependence : Steady-state or transient analysis
%                           .scheme : The time integration scheme
%                               .T0 : The start time of the simulation
%                             .TEnd : The end time of the simulation
%                          .nTSteps : The number of the time steps
%             intFEM : On the numerical integration (quadrature)
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
% 0. Read input
%
% 1. Solve the linear equation system
%
% 2. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_________________________________________________________\n');
    fprintf('#########################################################\n');
    fprintf('Static linear analysis for an isogeometric plate in plane\n');
    fprintf('stress action has been initiated\n');
    fprintf('_________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input
%FEM
% Number of nodes in the mesh
noNodesFEM = length(strMsh.nodes(:,1));

% Number of DOFs in the mesh
nDOFsFEM = 2*noNodesFEM;

% GLobal DOF numbering
DOFNumberingFEM = 1:nDOFsFEM;

% Assign dummy variables
uSaved = 'undefined';
uDot = 'undefined';
uDotSaved = 'undefined';
uDDotSaved = 'undefined';

% Title for the output file
title = 'geometrically linear steady-state plane stress analysis';

% Connectivity arrays for the DOFs into the resulting vectors
xDispFEM = 1:2:nDOFsFEM-1;
yDispFEM = 2:2:nDOFsFEM;

% Make directory to write out the results of the analysis
isExistent = exist(strcat('../../outputVTK/FEMPlateInMembraneActionAnalysis/',caseName),'dir');
if ~isExistent
    mkdir(strcat('../../outputVTK/FEMPlateInMembraneActionAnalysis/',caseName));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IGA
% Check input
mxi = length(Xi);
meta = length(Eta);
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Initialize dummy variables
dHatSaved = 'undefined';
dHatDotSaved = 'undefined';
dHatDDotSaved = 'undefined';
dHatDot = 'undefined';
nLinearAnalysis = 'undefined';
nLinearAnalysisFEM = 'undefined';
transientAnalysis = 'undefined';

% Assign the properties of the B-Spline patch
BSplinePatch.p = p;
BSplinePatch.q = q;
BSplinePatch.Xi = Xi;
BSplinePatch.Eta = Eta;
BSplinePatch.CP = CP;
BSplinePatch.isNURBS = isNURBS;

% Get number of Control Points
noCPs = length(CP(:,1,1))*length(CP(1,:,1));

% Get number of DOFs
noDOFs = 2*noCPs;

% Initialize output array
dHat = zeros(noDOFs+nDOFsFEM,1);

% Assign DoF numbering dof array assigns two DoF (x,y) to every CP 
% numbering follows CP: CP1->dof1,dof2 CP2->dof3,dof4

% Initialize the array of the degrees of freedome
DOFNumbering = zeros(nxi,neta,2);

% Initialize counter
k = 1;

% Loop over all the Control points
for cpj = 1:neta
    for cpi = 1:nxi
        DOFNumbering(cpi,cpj,1) = k;
        DOFNumbering(cpi,cpj,2) = k + 1;
        
        % Update counter
        k = k + 2;
    end
end

% Get a sequencial numbering of the unconstained DOFs into a vector
freeDOFs = zeros(length(Fl),1);
for i=1:length(Fl)
    freeDOFs(i,1) = i;
end
freeDOFsFEM = zeros(length(F),1);
for i=length(Fl)+1:length(Fl)+length(F)
    freeDOFs(i,1) = i;
end
freeDOFs(ismember(freeDOFs,homDOFs)) = [];

% The analysis is steady-state
t = 0;

%% 1. Solve the linear equation system
[dHat,~,~,FComplete,minElSize] = solve_IGA_FEM_LinearSystem(dHatSaved, ...
    dHatDotSaved,dHatDDotSaved,BSplinePatch,Fl,bodyForces,parameters,...
    dHat,dHatDot,@computeStiffMtxAndLoadVctIGAPlateInMembraneActionLinear,...
    DOFNumbering,freeDOFs,homDOFs,inhomDOFs,valuesInhomDOFs,transientAnalysis,...
    t,nLinearAnalysis,int,outMsg,...
    uSaved, ...
    uDotSaved,uDDotSaved,strMsh,analysis,F,bodyForcesFEM,parametersFEM,uDot, ...
    @computeStiffMtxAndLoadVctFEMPlateInMembraneActionCST,DOFNumberingFEM,...
    homDOFsFEM,inhomDOFsFEM,valuesInhomDOFsFEM,strDynamics,nLinearAnalysisFEM,intFEM,'');

%% 2. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Static linear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('______________Static Linear Analysis Ended_______________\n');
    fprintf('#########################################################\n\n\n');
end


