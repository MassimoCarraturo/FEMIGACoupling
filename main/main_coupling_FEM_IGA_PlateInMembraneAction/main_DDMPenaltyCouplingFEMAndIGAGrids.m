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
%% Script documentation
%
% Task : Benchmark of a quarter of a circle beam, fixed at the one edge and
%        subject to uniform shear pressure at the other end. Closed form 
%        solution exists in the polar coordinate system. The problem is
%        decomposed into two unsymmetric subdomains.
%
% Date : 26.11.2013
%
%% Preamble
clc;
clear;

%% Includes 
addpath('../../coupling_FEM_IGA/Solver/');
addpath('../../coupling_FEM_IGA/Solution/');
addpath('../../coupling_FEM_IGA/Graphics/');

% Add general math functions
addpath('../../generalMath/');

% Add all functions related to parsing
addpath('../../parsers/');

% Add all functions related to the low order basis functions
addpath('../../basisFunctions/');

% Equation system solvers
addpath('../../equationSystemSolvers/');

% Add all functions related to plate in membrane action analysis
addpath('../../FEMPlateInMembraneActionAnalysis/solvers/',...
        '../../FEMPlateInMembraneActionAnalysis/solutionMatricesAndVectors/',...
        '../../FEMPlateInMembraneActionAnalysis/loads/',...
        '../../FEMPlateInMembraneActionAnalysis/basisFunctions/',...
        '../../FEMPlateInMembraneActionAnalysis/graphics/',...
        '../../FEMPlateInMembraneActionAnalysis/output/',...
        '../../FEMPlateInMembraneActionAnalysis/postprocessing/');


% Add general auxiliary functions
addpath('../../auxiliary/');

% Add all functions related to the Computer-Aided Geometric Design (GACD) kernel
addpath('../../CAGDKernel/CAGDKernel_basisFunctions',...
        '../../CAGDKernel/CAGDKernel_geometryResolutionRefinement/',...
        '../../CAGDKernel/CAGDKernel_baseVectors/',...
        '../../CAGDKernel/CAGDKernel_graphics/',...
        '../../CAGDKernel/CAGDKernel_BsplineCurve/',...
        '../../CAGDKernel/CAGDKernel_BSplineSurface/');
    
% Add all functions related to the Isogeometric Plane Stress formulation
addpath('../../isogeometricPlateInMembraneActionAnalysis/stiffnessMatrices/',...
        '../../isogeometricPlateInMembraneActionAnalysis/loads/',...
        '../../isogeometricPlateInMembraneActionAnalysis/solvers/',...
        '../../isogeometricPlateInMembraneActionAnalysis/graphicsSinglePatch/',...
        '../../isogeometricPlateInMembraneActionAnalysis/graphicsMultipatches/',...
        '../../isogeometricPlateInMembraneActionAnalysis/postprocessing/',...
        '../../isogeometricPlateInMembraneActionAnalysis/auxiliary/',...
        '../../isogeometricPlateInMembraneActionAnalysis/errorComputation/',...
        '../../isogeometricPlateInMembraneActionAnalysis/lagrangeMultiplierDecomposition/',...
        '../../isogeometricPlateInMembraneActionAnalysis/penaltyDecomposition/',...
        '../../isogeometricPlateInMembraneActionAnalysis/perturbedLagrangeMultiplierDecomposition/',...
        '../../isogeometricPlateInMembraneActionAnalysis/partitionedGaussSeidelDecomposition/',...
        '../../isogeometricPlateInMembraneActionAnalysis/nitscheDecomposition/');
    
%% NURBS parameters

% Geometrical parameters
internalRadius1 = 4;
externalRadius1 = 5;
internalRadius2 = internalRadius1;
externalRadius2 = externalRadius1; 

% Polynomial degrees
p1 = 2;
q1 = 1;
p2 = 2;
q2 = 1;

% Knot vectors
Xi1 = [0 0 0 1 1 1];
Eta1 = [0 0 1 1];
Xi2 = [0 0 0 1 1 1];
Eta2 = [0 0 1 1];

% 1st patch :
% ___________

% This is modelled with the classical Finite Elements
% Define the path to the case
pathToCase = '../../inputGiD/FEM_IGA_TestCase/';
caseName = 'curvedBeamTipShear_new';

% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,IBC,analysis,parameters,nLinearAnalysis,strDynamics] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

% 2nd patch :
% ___________

% Control Point coordinates and weigthts

% x-coordinates
CP2(:,:,1) = [internalRadius2*cos(pi/4) externalRadius1*cos(pi/4)
             internalRadius2 externalRadius1
             internalRadius2 externalRadius1];
         
% y-coordinates
CP2(:,:,2) = [internalRadius2*cos(pi/4) externalRadius1*cos(pi/4);
             internalRadius2*(sqrt(2)-1) externalRadius1*(sqrt(2)-1)
             0 0];
         
% z-coordinates
CP2(:,:,3) = [0 0;
              0 0
              0 0];
          
% Control Point weights
w2 = cos(45*pi/180/2);
CP2(:,:,4) = [1 1;
              w2 w2
              1 1];
          
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS2 = 0;
nxi2 = length(CP2(:,1,1));
neta2 = length(CP2(1,:,1));
for i= 1:nxi2
    for j=1:neta2
        if CP2(i,j,4)~=1
            isNURBS2 = 1;
            break;
        end
    end
    if isNURBS2
        break;
    end
end

%% Material constants

% % 1st patch :
% % ___________
% 
% % Young's modulus
% parameters1.E = 1e5;
% 
% % Poisson ratio
% parameters1.nue = 0.333;
% 
% % Plate thickness
% parameters1.t = 1;
% 
% % 2nd patch :
% % ___________
% 
% % Young's modulus
% parameters2.E = 1e5;
% 
% % Poisson ratio
% parameters2.nue = 0.333;
% 
% % Plate thickness
% parameters2.t = 1;

%Penalty Term
Beta = 1e9;

%% GUI

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points

% 1st patch :
% ___________
% On the body forces
bodyForcesIGA = @computeConstantVecrticalBodyForceVct;

% On the numerical integration
% 'default', 'manual'
intLoad.type = 'default';
intDomain.type = 'default';
% Choice of the Gauss Point number
intLoad.nGP = 1;
intDomain.nGP = 1;

% Linear analysis
nLinearAnalysis = [];

% Initialize graphics index
graph.index = 1;


% 2nd patch :
% ___________
int2.type = 'default';
if strcmp(int2.type,'manual')
    int2.xiNGP = 6;
    int2.etaNGP = 6;
    int2.xiNGPForLoad = 6;
    int2.xetaNGPForLoad = 6;
end

% Interface integration :
% _______________________
intC.type = 'default';
intC.method = 'lagrangeMultipliers';
if strcmp(intC.type,'manual')
    if strcmp(intC.method,'lagrangeMultipliers')
        intC.nGP1 = 6;
        intC.nGP2 = 6;
    else
        intC.nGP = 6;
    end
end

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement', 'strain', 'stress'
graph.resultant = 'displacement';

% Component of the resultant to plot
% .component: 'x', 'y', '2norm','xy', '1Principal', '2Principal'
graph.component = '2norm';

% error computation for this benchmark
% .resultant: 'displacement', 'strain', 'stress'
error.resultant = 'stress';
% .component: 'x', 'y', 'tensor','xy'
error.component = 'tensor';
% .xiNGP .eta.NGP
error.xiNGP = 10;
error.etaNGP = 10;

%% Refinements

% polynomial order
a = 1;
b = 1;

% 1st patch :
% ___________

% 2nd patch :
% ___________

% Degree elevation 
tp2 = b;  tq2 = b;
[Xi2,Eta2,CP2,p2,q2] = degreeElevateBSplineSurface(p2,q2,Xi2,Eta2,CP2,tp2,tq2,'outputEnabled');

% Knot insertion
n2 = ceil(7*2);
[Xi2,Eta2,CP2] = knotRefineUniformlyBSplineSurface(p2,Xi2,q2,Eta2,CP2,n2,n2,'outputEnabled');

lng = length(Xi2);
XiB = Xi2(1);

%% Dirichlet and Neumann boundary conditions

% 1st patch :
% ___________
homDOFs=[];
inhomDOFs=[];
valuesInhomDOFs=[];

% 2nd patch :
% ___________

% Supports
rb2 = [];
% xisup = [1 1];   etasup = [0 1];   directionSupport = 2;
% rb2 = findDofs2D(rb2,usup,vsup,dirs,CP2);

% Load
Fl2 = [];
fx = -1;
xib = 1;   etab = [0 1];   directionForce=1;
Fl2 = computeLoadVctLineIGAPlateInMembraneAction(Fl2,xib,etab,p2,q2,Xi2,Eta2,CP2,isNURBS2,fx,directionForce,int2,'outputEnabled');
% xib = [0 0.5];   etab = 1;   directionForce = 1;
% fl2 = computeLoadVctLineIGAPlateInMembraneAction(Fl2,xib,etab,p2,q2,Xi2,Eta2,CP2,isNURBS2,fx,directionForce,int2,'outputEnabled');

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Compute the load vector FEM
bodyForcesFEM = @computeConstantVecrticalBodyForceVct;
t = 0;
F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,intLoad,'outputEnabled');
%% Find the degrees of freedom over the coupled interfaces

% 1st patch :
% ___________


% 2nd patch :
% ___________

cb2 = [];
xicoup2 = [0 0];   
etacoup2 = [0 1];

%% Projection

% P = strMsh.nodes(2,:)';
% 
% [xiP,etaP,Projected,flag,noIterations] = ...
%                         computeNearestPointProjectionOnBSplineSurface...
%                         (P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,newtonRaphson)


%% Solve the system
[dHat,FComplete,minElSize,K] = solve_IGA_FEM_PlateInMembraneActionLinear(Beta,XiB,p2,Xi2,q2,Eta2,CP2,isNURBS2,parameters,Fl2,bodyForcesIGA,homDOFs,...
    inhomDOFs,valuesInhomDOFs,int2,'outputEnabled',strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,IBC,F,bodyForcesFEM,...
    analysis,parameters,nLinearAnalysis,strDynamics,intDomain,caseName,pathToOutput,'outputEnabled');

%% Postprocessing
noDOFsFEM = length(F);
homDOFsAndreas = [homDBC homDOFs+noDOFsFEM+1];

% Energy norm
Energynorm=dHat'*K*dHat;
Energynorm=sqrt(0.5*Energynorm);

Kred = K;
Kred(homDOFsAndreas,:) = [];
Kred(:,homDOFsAndreas) = [];
 maxEigenValue = max(eigs(Kred));
 minEigenValue = min(eigs(Kred,6,'sm'));
 conditionNumber = maxEigenValue/minEigenValue;

% Plotting
graph.visualization.geometry = 'reference_and_current';
noDOFsFEM = length(strMsh.nodes(:,1))*2;
graph.index = plot_FEM_IGA_currentConfigurationAndResultants(p2,q2,Xi2,Eta2,CP2,isNURBS2,homDOFs,parameters,Fl2,dHat(noDOFsFEM+1:length(dHat)),graph,'outputEnabled',...
    strMsh,homDBC,dHat(1:noDOFsFEM),parameters,analysis);




%% End