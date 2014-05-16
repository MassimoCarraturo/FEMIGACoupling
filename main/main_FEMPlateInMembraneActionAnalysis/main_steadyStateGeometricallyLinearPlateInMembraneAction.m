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
% Task : Plane stress analysis for a rectangular plate subject to uniform
%        pressure on its top edge
%
% Date : 19.02.2014
%
%% Preamble

% Clear memory
clear;

% Clear the command window
clc;

% Close all generated windows
close all;

%% Includes

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

%% Parse data from GiD input file

% Define the path to the case
pathToCase = '../../inputGiD/FEMPlateInMembraneActionAnalysis/';
% caseName = 'cantileverBeamPlaneStress';
% caseName = 'PlateWithAHolePlaneStress';
% caseName = 'PlateWithMultipleHolesPlaneStress';
% caseName = 'InfinitePlateWithAHolePlaneStress';
caseName = 'curvedPlateTipShearPlaneStress';

% Parse the data from the GiD input file
[strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,nLinearAnalysis,strDynamics] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,'outputEnabled');

%% GUI

% On the body forces
bodyForces = @computeConstantVecrticalBodyForceVct;

% On the numerical integration
% 'default', 'manual'
intLoad.type = 'manual';
intDomain.type = 'manual';
% Choice of the Gauss Point number
intLoad.nGP = 1;
intDomain.nGP = 1;

% Linear analysis
nLinearAnalysis = [];

% Initialize graphics index
graph.index = 1;

%% Output data to a VTK format
pathToOutput = '../../outputVTK/FEMPlateInMembraneActionAnalysis/';

%% Compute the load vector
t = 0;
F = computeLoadVctFEMPlateInMembraneAction(strMsh,analysis,NBC,t,intLoad,'outputEnabled');

%% Visualization of the configuration
graph.index = plot_referenceConfigurationFEMPlateInMembraneAction(strMsh,analysis,F,homDBC,graph,'outputEnabled');

%% Solve the plate in membrane action problem
[dHat,FComplete,minElSize] = ...
    solve_FEMPlateInMembraneAction(strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,F,bodyForces,analysis,parameters,nLinearAnalysis,strDynamics,intDomain,caseName,pathToOutput,'outputEnabled');

%% Postprocessing
graph.visualization.geometry = 'reference_and_current';
resultant = 'stress';
component = 'y';
graph.index = plot_currentConfigurationAndResultants(strMsh,homDBC,dHat,parameters,analysis,resultant,component,graph);

%% END OF THE SCRIPT