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
function [dHat1,dHat2,lambdaHat,F1,F2,rankDeficiency,condK] = solve_DDMLagrangeMultiplierIGAPlateInMembraneAction(patch1,patch2,lagrangeMultiplierField,intC,outMsg)
%% Function documentation
%
% Returns the displacement and the Lagange Mutlipliers fields as well as 
% the complete force vectors corresponding to the Lagrange Multiplier 
% decomposition method for its application to the isogeometric plate in 
% membrane action analysis decomposed into two subdomains.
%
%                  Input :
%          patch1,patch2 : Structures containing all geometrical 
%                          information (NURBS parameters), technical 
%                          information (parameters) and the coupling 
%                          information (structure cb)
% lagrangeMultiplierField : Structure containing all information on 
%                           polynomial order and number of elements for the 
%                           Lagrange multipliers field on the coupling 
%                           surface
%                 Fl1,Fl2 : The load vectors corresponding to each patch
%                           seperately
%                    intC : Structure containing information on the
%                           numerical integration over the coupling surface 
%                  outMsg : Whether or not to output message on refinement 
%                           progress 'outputEnabled' : enables output 
%                           information
%   
%                  Output :
%             dHat1,dHat2 : The complete displacement vectors for the two
%                           patches
%               lambdaHat : The discrete Lagrange multipliers vector
%                   F1,F2 : The complete force vectors for both patches
%          rankDeficiency : The rank deficiency of the saddle point system
%                   condK : Condition number of the reduced equation system
%
% Function Layout :
%
% 1. Compute the stiffness matrices for the two patches
%
% 2. Determine the direction of the coupling interface
%
% 3. Compute the coupling matrices for the two patches
%
% 4. Assemble to the global stiffness matrix
%
% 5. Compute the generalized load vector
%
% 6. Reduce the system according to the given constraints
%
% 7. Miscellaneous operations
%
% 8. Solve for the global displacement vector
%
% 9. Re-assign the output arrays
%
% 10. Compute the complete force vectors
% 
% 11. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________\n');
    fprintf('#############################################################\n');
    fprintf('Static linear analysis for the decomposed into two subdomains \n');
    fprintf('isogeometric plate in plane stress action using the Lagrange \n');
    fprintf('Multiplier method has been initiated \n\n');
    fprintf('Discretization for patch 1 :\n');
    fprintf('--------------------------- \n');
    fprintf('Polynomial degree in xi-direction = %d \n',patch1.p);
    fprintf('Polynomial degree in eta-direction = %d \n',patch1.q);
    fprintf('Number of elements = %d \n\n',(length(patch1.Xi) - 2*(patch1.p + 1)+1)*...
        (length(patch1.Eta) - 2*(patch1.q + 1)+1));
    fprintf('Discretization for patch 2 :\n');
    fprintf('--------------------------- \n');
    fprintf('Polynomial degree in xi-direction = %d \n',patch2.p);
    fprintf('Polynomial degree in eta-direction = %d \n',patch2.q);
    fprintf('Number of elements = %d \n \n',(length(patch2.Xi) - 2*(patch2.p + 1)+1)*...
        (length(patch2.Eta) - 2*(patch2.q + 1)+1));
    fprintf('Discretization for the Lagrange Multiplier Field :\n');
    fprintf('-------------------------------------------------- \n');
    fprintf('Polynomial degree = %d \n',lagrangeMultiplierField.p);
    fprintf('Number of elements = %d \n',length(lagrangeMultiplierField.Xi) - 2*(lagrangeMultiplierField.p + 1) + 1);
    fprintf('______________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% 1st patch :
% ___________

p1 = patch1.p;
q1 = patch1.q;
Xi1 = patch1.Xi;
Eta1 = patch1.Eta;
CP1 = patch1.CP;
isNURBS1 = patch1.isNURBS;
parameters1 = patch1.parameters;
rb1 = patch1.rb;
Fl1 = patch1.Fl;
xicoup1 = patch1.xicoup;
etacoup1 = patch1.etacoup;
int1 = patch1.int;

% 2nd patch :
% ___________

p2 = patch2.p;
q2 = patch2.q;
Xi2 = patch2.Xi;
Eta2 = patch2.Eta;
CP2 = patch2.CP;
isNURBS2 = patch2.isNURBS;
parameters2 = patch2.parameters;
rb2 = patch2.rb;
Fl2 = patch2.Fl;
xicoup2 = patch2.xicoup;
etacoup2 = patch2.etacoup;
int2 = patch2.int;

%% 1. Compute the stiffness matrices for the two patches
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the stiffness matrices for both patches\n');
end

% 1st patch :
% ___________

K1 = computeStiffMtxIGAPlateInMembraneActionLinear(p1,Xi1,q1,Eta1,CP1,isNURBS1,parameters1,int1,'');
[s1,~] = size(K1);
isMaster = 1;

% 2nd patch :
% ___________

K2 = computeStiffMtxIGAPlateInMembraneActionLinear(p2,Xi2,q2,Eta2,CP2,isNURBS2,parameters2,int2,'');
[s2,~] = size(K2);
isSlave = 0;

%% 2. Determine the direction of the coupling interface
haveTheSameOrientation = findSubdomainInterfaceOrientation(p1,Xi1,q1,Eta1,CP1,isNURBS1,xicoup1,etacoup1,p2,Xi2,q2,Eta2,CP2,isNURBS2,xicoup2,etacoup2);

%% 3. Compute the coupling matrices for the two patches
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the Lagrange Multiplier coupling matrices for both patches\n');
end

% 1st patch :
% ___________

C1 = computeLagrangeMultiplierCouplingMatrixIGAPlateInMembraneAction(patch1,lagrangeMultiplierField,isMaster,haveTheSameOrientation,intC);
[~,sm12] = size(C1);

% 2nd patch :
% ___________

C2 = computeLagrangeMultiplierCouplingMatrixIGAPlateInMembraneAction(patch2,lagrangeMultiplierField,isSlave,haveTheSameOrientation,intC);

%% 4. Assemble to the global stiffness matrix
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Assembling to the global Lagrange Multiplier matrix of the coupled system \n');
end

% Form the null matrices
nullMatrix12 = zeros(s1,s2);
nullMatrix33 = zeros(sm12,sm12);

% Global system matrix for the coupled system
K = [K1            nullMatrix12 C1
     nullMatrix12' K2           C2 
     C1'           C2'          nullMatrix33];

%% 5. Compute the generalized load vector

% Form the null matrice
nullMatrix3 = zeros(sm12,1);

% Generalized load vector
F = [Fl1
     Fl2 
     nullMatrix3];

%% 6. Reduce the system according to the given constraints
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Reducing the system according to the given constraints\n');
end

% Rearrange the support vectors with respect to the sequential numbering of
% the full system
for i=1:length(rb2)
    rb2(i) = s1 + rb2(i);
end

% Merge the supports vectors into one and delete double entries if any
rb = mergesorted(rb1,rb2);
rb = unique(rb);

% Reduce the stiffness matrix and the load vector
Kred = K;
Fred = F;
Kred(:,rb) = [];
Kred(rb,:) = [];
Fred(rb) = [];

%% 7. Miscellaneous operations
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing rank and conditioning of the reduced system matrix \n');
    
    % Size of the reduced system matrix
    [sizeKred,~] = size(Kred);

    % Rank of the reduced system matrix
    rankKred = rank(Kred);

    % Conditioning of the reduced matrix
    condK = cond(Kred);

    % Message on the rank deficiency of the reduced system matrix
    if sizeKred==rankKred
        rankDeficiency = sizeKred-rankKred;
        fprintf('>> Inf-Sup conditions are satisfied, the reduced system has full rank\n');
    else
        rankDeficiency = sizeKred-rankKred;
        fprintf('>> Warning : Inf-Sup conditions are violated, the reduced system has rank deficiency equal to %d \n',rankDeficiency);
    end
else
    rankDeficiency = 'undefined';
    condK = 'undefined';
end

%% 8. Solve for the global displacement vector
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Solving the linear system of %d equations\n',length(Fred));
end
Kred = sparse(Kred);
dHatred = Kred\Fred;

%% 9. Re-assign the output arrays
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Re-assembling to the complete displacement vector\n');
end

% Assemble complete displacement vector
% Dimension of the complete displacement vector
dHat = zeros(length(F),1);

% Get a sequencial numbering of the unconstained DOFs into a vector
unconstrainedDOFs = zeros(length(F),1);
for i=1:length(F)
    unconstrainedDOFs(i,1) = i;
end
unconstrainedDOFs(ismember(unconstrainedDOFs,rb)) = [];

% Fill up the complete displacement field
dHat(unconstrainedDOFs) = dHatred;

% Assign the values to each displacement vector respectively
dHat1 = dHat(1:s1);
dHat2 = dHat(s1+1:s1+s2);
lambdaHat = dHat(s1+s2+1:length(dHat));

%% 10. Compute the complete force vectors
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the complete force vectors\n\n');
end
F1 = K1*dHat1;
F2 = K2*dHat2;

%% 11. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Static linear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('_________________Static Linear Analysis Ended_________________\n');
    fprintf('##############################################################\n\n\n');
end

end
