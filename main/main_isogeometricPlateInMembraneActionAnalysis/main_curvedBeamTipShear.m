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
% Task : The benchmark example is a curved beam (modelled as a thin plate) 
%        which is subject into tip shear force. For this example there is 
%        closed form solution in terms of the stress field
%
% Date : 15.11.2013
%
%% Preamble
clc;
clear;

%% Includes

% Add general math functions
addpath('../../generalMath/');

% Add general auxiliary functions
addpath('../../auxiliary/');

% Add linear solvers
addpath('../../equationSystemSolvers/');

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
        '../../isogeometricPlateInMembraneActionAnalysis/postprocessing/',...
        '../../isogeometricPlateInMembraneActionAnalysis/auxiliary/',...
        '../../isogeometricPlateInMembraneActionAnalysis/errorComputation/');
    
%% NURBS parameters

% Geometrical parameters
internalRadius = 4;
externalRadius = 5;

% Polynomial degrees
p = 2;
q = 1;

% Knot vectors
Xi = [0 0 0 1 1 1];
Eta = [0 0 1 1];

% Control Point coordinates

% x-coordinates
CP(:,:,1) = [0 0
             internalRadius externalRadius;
             internalRadius externalRadius];
         
% y-coordinates
CP(:,:,2) = [internalRadius externalRadius
             internalRadius externalRadius
             0 0];
         
% z-coordinates
CP(:,:,3) = [0 0;
             0 0
             0 0];

% Weights
w = cos(45*pi/180);
CP(:,:,4) = [1 1;
              w w
              1 1];
          
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS = 0;
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));
for i= 1:nxi
    for j=1:neta
        if CP(i,j,4)~=1
            isNURBS = 1;
            break;
        end
    end
    if isNURBS
        break;
    end
end

%% Material constants 

% Young's modulus
parameters.E = 1e5;

% Poisson ratio
parameters.nue = 0.0;

% Thickness of the plate
parameters.t = 1;

%% GUI

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'manual' : manual choice of the number of Gauss points
int.type = 'default';
if strcmp(int.type,'manual')
    int.xiNGP = 6;
    int.etaNGP = 6;
    int.xiNGPForLoad = 6;
    int.xetaNGPForLoad = 6;
end

% On the graphics
graph.index = 1;

% On the postprocessing:
% .postprocConfig : 'reference','current','referenceCurrent'
graph.postprocConfig = 'referenceCurrent';

% Plot strain or stress field
% .resultant: 'displacement', 'strain', 'stress'
graph.resultant = 'stress';

% Component of the resultant to plot
% .component: 'x', 'y', '2norm','xy', '1Principal', '2Principal'
graph.component = 'xy';

% error computation for this benchmark
% .resultant: 'displacement', 'strain', 'stress'
error.resultant = 'stress';
% .component: 'x', 'y', 'tensor','xy'
error.component = 'tensor';
% .xiNGP .eta.NGP
error.xiNGP = 16;
error.etaNGP = 16;

% Function handle to the computation of the body force vector
bodyForces = @computeConstantVecrticalBodyForceVct;

%% Refinement 

% Degree elevation of the surface
a = 0;
tp = a;  tq = a;
[Xi,Eta,CP,p,q] = degreeElevateBSplineSurface(p,q,Xi,Eta,CP,tp,tq,'outputEnabled');

% Knot insertion on the surface
a = 2;
refXi = ceil(25*a);
refEta = ceil(7*a);
[Xi,Eta,CP] = knotRefineUniformlyBSplineSurface(p,Xi,q,Eta,CP,refXi,refEta,'outputEnabled');

%% Dirichlet and Neumann boundary conditions 

% homogeneous Dirichlet boundary conditions
homDOFs = [];
xiSup = [0 0];   etaSup = [0 1];   dirSupp = 1;
homDOFs = findDofs2D(homDOFs,xiSup,etaSup,dirSupp,CP);
xiSup = [0 0];   etaSup = [0 0];   dirSupp = 2;
homDOFs = findDofs2D(homDOFs,xiSup,etaSup,dirSupp,CP);

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs = [];
valuesInhomDOFs = [];

% load (Neuman boundary conditions)
FAmp = -1;  
Fl = [];
xib = 1;   etab = [0 1];   dirFor = 1;
Fl = computeLoadVctLineIGAPlateInMembraneAction(Fl,xib,etab,p,q,Xi,Eta,CP,isNURBS,FAmp,dirFor,int,'outputEnabled');

%% Plot the initial configuration
figure(graph.index)
plot_referenceConfigurationIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,homDOFs,Fl,'outputEnabled');
title('Reference configuration for an isogeometric plate in membrane action');
graph.index = graph.index + 1;

%% Solve the system
[dHat,~,K] = solve_IGAPlateInMembraneActionLinear(p,Xi,q,Eta,CP,isNURBS,parameters,Fl,bodyForces,homDOFs,...
    inhomDOFs,valuesInhomDOFs,int,'outputEnabled');

%% Postprocessing 
% Energy norm
Energynorm=dHat'*K*dHat;
Energynorm=sqrt(0.5*Energynorm);

% Plot the current configuration and the selected resultant
%graph.index = plot_postprocIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,homDOFs,parameters,Fl,dHat,graph,'outputEnabled');

% Compute the relative error in the L2-norm for the stress resultant
[errorCurvedBeamTipShear,minElArea] = computeRelErrorL2CurvedBeamTipShearIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,parameters,internalRadius,externalRadius,abs(FAmp),dHat,error,'outputEnabled');

%% End of the scrpit