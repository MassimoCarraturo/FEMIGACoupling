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
    
% This is modelled with the classical Finite Elements
% Define the path to the case
pathToCase = '../../inputGiD/FEM_IGA_TestCase/';
caseName = 'quarterPlate2';

% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,IBC,analysis,parameters,nLinearAnalysis,strDynamics] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

    %% NURBS parameters

% Geometrical parameters of the problem
R = 1;
T = 10;

% Polynomial degrees
p1 = 1;
q1 = 2;
p2 = 2;
q2 = 1;

% Knot vectors
Xi1 = [0 0 1 1];
Eta1 = [0 0 0 1 1 1];

Xi2 = [0 0 0 1 1 1];
Eta2 = [0 0 1 1];

% 1st patch :
% ___________

% x-coordinates
CP1(:,:,1) = [-4 -4 -4;
              -1 -1 -sqrt(2)/2];

% y-coordinates
CP1(:,:,2) = [0 2                   4
              0 (2-sqrt(2))/sqrt(2) sqrt(2)/2];

% z-coordinates
CP1(:,:,3) = [1 1 1 
              1 1 1];

% Control Point weights
w1 =  cos(22.5*pi/180);
CP1(:,:,4) = [1 1 1; 
              1 w1 1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = 0;
nxi1 = length(CP1(:,1,1));
neta1 = length(CP1(1,:,1));
for i= 1:nxi1
    for j=1:neta1
        if CP1(i,j,4)~=1
            isNURBS1 = 1;
            break;
        end
    end
    if isNURBS1
        break;
    end
end


%% Material constants

% 1st patch :
% ___________

% On the body forces
bodyForcesIGA = @computeConstantVecrticalBodyForceVct;

% Young's modulus
parameters1.E = 1e5;

% Poisson ratio
parameters1.nue = 0.0;

% Plate thickness
parameters1.t = 1;


%% GUI

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points

% 1st patch :
% ___________
int1.type = 'default';
if strcmp(int1.type,'manual')
    int1.xiNGP = 12;
    int1.etaNGP = 12;
    int1.xiNGPForLoad = 12;
    int1.xetaNGPForLoad = 12;
end

% 2nd patch :
% ___________
int2.type = 'default';
if strcmp(int2.type,'manual')
    int2.xiNGP = 12;
    int2.etaNGP = 12;
    int2.xiNGPForLoad = 12;
    int2.xetaNGPForLoad = 12;
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
errorDisplacements.resultant = 'displacement';
errorStress.resultant = 'stress';
% .component: 'x', 'y', 'tensor','xy'
errorDisplacements.component = 'tensor'; 
errorStress.component = 'tensor';
% .xiNGP .eta.NGP
errorDisplacements.xiNGP = 10;
errorDisplacements.etaNGP = 10;
errorStress.xiNGP = 10;
errorStress.etaNGP = 10;

%% Refinements

% polynomial order
a = 1;
b = 0;

% 1st patch :
% ___________

% Degree elevation 
tp1 = a;  tq1 = a;
[Xi1,Eta1,CP1,p1,q1] = degreeElevateBSplineSurface(p1,q1,Xi1,Eta1,CP1,tp1,tq1,'outputEnabled');

% Knot insertion
scale = 1;
n1 = ceil(3*scale);
[Xi1,Eta1,CP1] = knotRefineUniformlyBSplineSurface(p1,Xi1,q1,Eta1,CP1,n1,n1,'outputEnabled');

lng = length(Xi1);
XiB = Xi1(lng);

%% Dirichlet and Neumann boundary conditions

% 1st patch :
% ___________

% Supports
rb1 = [];
xisup = [0 1];   etasup = [0 0];   directionSupport = 2;
rb1 = findDofs2D(rb1,xisup,etasup,directionSupport,CP1);

% Load
Fl1 = [];
fx1 = @computeSigmaXXForAnalyticalSolutionToInfinitePlateWithHole;
fy1 = @computeSigmaXYForAnalyticalSolutionToInfinitePlateWithHole;
xib1 = 0;   etab1 = [0 1];
Fl1 = computeLoadVctLineIGAPlateInMembraneAction(Fl1,xib1,etab1,p1,q1,Xi1,Eta1,CP1,isNURBS1,fx1,1,int1,'');
Fl1 = computeLoadVctLineIGAPlateInMembraneAction(Fl1,xib1,etab1,p1,q1,Xi1,Eta1,CP1,isNURBS1,fy1,2,int1,'');
Fl1 = - Fl1;


%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Compute the load vector FEM
bodyForcesFEM = @computeConstantVecrticalBodyForceVct;
t = 0;
% On the numerical integration
% 'default', 'manual'
intLoad.type = 'default';
intDomain.type = 'default';

F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,intLoad,'outputEnabled');
%% Find the degrees of freedom over the coupled interfaces

% 1st patch :
% ___________
homDOFs=[];
inhomDOFs=[];
valuesInhomDOFs=[];

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

%Penalty Term
Beta = 1e9;

%% Solve the system
[dHat,FComplete,minElSize,K] = solve_IGA_FEM_PlateInMembraneActionLinear(Beta,XiB,p1,Xi1,q1,Eta1,CP1,isNURBS1,parameters,Fl1,bodyForcesIGA,homDOFs,...
    inhomDOFs,valuesInhomDOFs,int2,'outputEnabled',strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,IBC,F,bodyForcesFEM,...
    analysis,parameters,nLinearAnalysis,strDynamics,intDomain,caseName,pathToOutput,'outputEnabled');

%% Postprocessing
noDOFsFEM = length(F);
homDOFsAndreas = [homDBC homDOFs+noDOFsFEM+1];

% Energy norm
Energynorm=dHat'*K*dHat;
Energynorm=sqrt(0.5*Energynorm);
% 
% Kred = K;
% Kred(homDOFsAndreas,:) = [];
% Kred(:,homDOFsAndreas) = [];
%  maxEigenValue = max(eigs(Kred));
%  minEigenValue = min(eigs(Kred,6,'sm'));
%  conditionNumber = maxEigenValue/minEigenValue;

% Plotting
graph.visualization.geometry = 'reference_and_current';
noDOFsFEM = length(strMsh.nodes(:,1))*2;
graph.index = plot_FEM_IGA_currentConfigurationAndResultants(p1,q1,Xi1,Eta1,CP1,isNURBS1,homDOFs,parameters,Fl1,dHat(noDOFsFEM+1:length(dHat)),graph,'outputEnabled',...
    strMsh,homDBC,dHat(1:noDOFsFEM),parameters,analysis);




%% End