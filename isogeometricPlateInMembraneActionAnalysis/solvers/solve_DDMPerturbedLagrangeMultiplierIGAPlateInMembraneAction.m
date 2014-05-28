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
function [dHat1,dHat2,lambdaHat,F1,F2,rDLambda,cNLambda,rDCLamdaC,cNCLamdaC] = solve_DDMPerturbedLagrangeMultiplierIGAPlateInMembraneAction(patch1,...
    patch2,lagrangeMultiplierField,alpha,epsilon,fixedPointTolerance,relaxationFactor,maxNoIterations,intC,outMsg)
%% Function documentation
%
% Returns the displacement and the Lagange Mutlipliers fields as well as 
% the complete force vectors corresponding to the perturbed Lagrange 
% Multiplier decomposition method for its application to the isogeometric 
% plate in membrane action analysis decomposed into two subdomains. Within
% the perturbed Lagrange Multiplier approach a small quadratic term is
% added into the functional such that the positive definiteness of the
% system matrices is ensured, plus the same quadratic term is added on the
% right-hand side of the equations such that fixed point iterations ensure
% the fulfilment of the Dirichlet compatibility conditions on the 
% interface. Additionally, penalty terms are also apparent in the
% formulation which can be neglected in case the penalty paramater is
% chosen to be zero.
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
%                   alpha : The penalty parameter
%                 epsilon : The perturbation factor
%     fixedPointTolerance : Convergence tolerance for the fixed point
%                           iterations
%        relaxationFactor : Relaxation factor for the fixed point
%                           iterations, relaxationFactor \in (0,1]
%         maxNoIterations : Maximum number of fixed point iterations
%                    intC : Structure containing information on the
%                           numerical integration over the coupling surface 
%                  outMsg : Whether or not to output message on refinement 
%                           progress 'outputEnabled' : enables output 
%                           information
%
%                  Output :
%             dHat1,dHat2 : The Control Point displacement field for both
%                           patches
%               lambdaHat : The Control Point displacement field for the
%                           Lagrange Multiplier field
%                   F1,F2 : The complete force vectors for both patches
%      rDLambda,rDCLamdaC : The rank deficiency for both matrices Lambda
%                           and [C1' C2']*Lamda*[C1; C2]
%      cNLambda,cNCLamdaC : The condition number for both matrices Lambda
%                           and [C1' C2']*Lamda*[C1; C2]
%
% Function Layout :
%
% 0. Read input
%
% 1. Compute the stiffness matrices for the two patches
%
% 2. Determine the direction of the coupling interface
%
% 3. Compute the coupling matrices for the two patches
%
% 4. Compute the penalty contributions to the stiffness matrices
%
% 5. Compute the penalty contributions to the coupling matrices
%
% 6. Compute the perturbation matrix
%
% 7. Assemble to the global matrix
%
% 8. Compute the generalized load vector
%
% 9. Reduce the system according to the constraints
%
% 10. Miscellaneous operations
%
% 11. Loop over all the fixed point iterations
% ->
%     11i. Update the force vector
%
%    11ii. Solve the structural problem
%
%   11iii. Solve the perturbed interface coupling equation
%
%    11iv. Compute the updated Lagrange Multiplier field
%
%     11v. Compute the error in the 2-norm for the increment in the Lagrange Multipliers field
%
%    11vi. Print message on the residual error
%
%   11vii. Update iteration counter
% <-    
% 12. Print message on the convergence of the fixed point iterations
%
% 13. Re-assign the output arrays
%
% 14. Compute the complete force vectors for both patches
%
% 15. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________________________\n');
    fprintf('###################################################################################\n');
    fprintf('Static linear analysis for the decomposed into two subdomains isogeometric plate in \n');
    fprintf('plane stress action using the perturbed Lagrange Multiplier method has been \n');
    fprintf('initiated \n\n');
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
    fprintf('Number of elements = %d \n \n',length(lagrangeMultiplierField.Xi) - 2*(lagrangeMultiplierField.p + 1) + 1);
    fprintf('Perturbation and augmentation parameters :\n');
    fprintf('------------------------------------------\n');
    fprintf('Penalty parameter alpha = %d\n',alpha);
    fprintf('Perturbation factor epsilon = %d\n\n',epsilon);
    fprintf('Fixed Point iteration parameters :\n');
    fprintf('----------------------------------\n');
    fprintf('Fixed point iteration tolerance = %d\n',fixedPointTolerance);
    fprintf('Under-relaxation factor = %d\n',relaxationFactor);
    fprintf('Maximum number of fixed point iterations = %d\n',maxNoIterations);
    fprintf('_________________________________________________________________________________\n\n');

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

% Get an initial value for the error
errordLambda = 2*fixedPointTolerance;

% Initialize counter for the fixed point iterations
counterFixedPointIterations = 1;

% Initialize string
if strcmp(outMsg,'outputEnabled')
    reverseStrFixedPointIterations = '';
end

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

% For patch 2:
C2 = computeLagrangeMultiplierCouplingMatrixIGAPlateInMembraneAction(patch2,lagrangeMultiplierField,isSlave,haveTheSameOrientation,intC);

% Form the global coupling matrix
C = [C1
     C2];
 
%% 4. Compute the penalty contributions to the stiffness matrices
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the penalty contributions to the stiffness matrices for both patches\n');
end

% 1st patch :
% ___________

if alpha~=0
    KPenalty1 = computePenaltyStiffnessMatrixIGAPlateInMembraneAction(patch1,alpha,isMaster,haveTheSameOrientation,intC);
else
    KPenalty1 = zeros(2*length(CP1(:,1,1))*length(CP1(1,:,1)),2*length(CP1(:,1,1))*length(CP1(1,:,1)));
end

% 2nd patch :
% ___________

if alpha~=0
    KPenalty2 = computePenaltyStiffnessMatrixIGAPlateInMembraneAction(patch2,alpha,isSlave,haveTheSameOrientation,intC);
else
    KPenalty2 = zeros(2*length(CP2(:,1,1))*length(CP2(1,:,1)),2*length(CP2(:,1,1))*length(CP2(1,:,1)));
end

%% 5. Compute the penalty contributions to the coupling matrices
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the penalty contributions to the coupling matrices for both patches\n');
end

% 1st patch :
% ___________

if alpha~=0
    CPenalty1 = computePenaltyCouplingMatrixIGAPlateInMembraneAction(patch1,patch2,alpha,isMaster,haveTheSameOrientation,intC);
else
    CPenalty1 = zeros(2*length(CP1(:,1,1))*length(CP1(1,:,1)),2*length(CP2(:,1,1))*length(CP2(1,:,1)));
end

% 2nd patch :
% ___________

% No further computation is needed for the penalty contribution to the
% coupling matrix on the slave side cause the transpose of CPenalty1 can be
% used instead due to symmetry

%% 6. Compute the perturbation matrix
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the perturbation matrix\n');
end
Lambda = computePerturbationMatrixIGAPlateInMembraneAction(lagrangeMultiplierField,patch1,epsilon,intC);

% In what concerns the discrete vector of unknowns for the Lagrange 
% Multipliers field, the idea is taken by Felippa's script; When the 
% iterations converge the interface Dirichlet compatibility condition 
% should be fulfiled up to the given tolerance for the fixed point 
% iterations

% Choose a predictor to be equal to zero
lambdaHat = zeros(length(Lambda),1);

%% 7. Assemble to the global matrix
if strcmp(outMsg,'outputEnabled')
    if alpha~=0
        fprintf('>> Assembling to the master penalty-coupling matrix\n');
    else
        fprintf('>> Assembling to the master stiffness matrix\n');
    end
end
K = [K1 + KPenalty1  CPenalty1
     CPenalty1'      K2 + KPenalty2];

%% 8. Compute the generalized load vector
F = [Fl1
     Fl2];

%% 9. Reduce the system according to the constraints
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

% Compute the product perturbation matrix
CLambdaC = C*(Lambda\C');

% Reduce the system matrices with respect to the constraints
Fred = F;
Kred = K;

% Reduce the stiffness matrix
Kred(:,rb) = [];
Kred(rb,:) = [];

% Reduce the product perturbation matrix
CLambdaC(:,rb) = [];
CLambdaC(rb,:) = [];

% Reduce the coupling matrix
C(rb,:) = [];

% Reduce the load vector
Fred(rb) = [];

%% 10. Miscellaneous operations
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing rank and conditioning of the reduced system matrices \n');
    
    % Compute the rank deficiency of matrix Lambda
    [sLambda,~] = size(Lambda);
    rDLambda = sLambda - rank(Lambda);
    if rDLambda~=0
        fprintf('>> Perturbation matrix Lambda has rank deficiency equal to %d \n',rDLambda);
    else
        fprintf('>> Perturbation matrix Lambda has full rank \n');
    end
    
    % Compute the condition number of matrix Lambda
    cNLambda = cond(Lambda);
    fprintf('>> Perturbation matrix has condition number equal to %d \n',cNLambda);
    
    % Compute the rank deficiency of matrix Kred-MLM
    [sKCLamdaC,~] = size(Kred-CLambdaC);
    rDCLamdaC = sKCLamdaC - rank(Kred-CLambdaC);
        if rDLambda~=0
        fprintf('>> Matrix Kred-MLM has rank deficiency equal to %d \n',rDCLamdaC);
    else
        fprintf('>> Matrix Kred-MLM has full rank\n');
        end
    
    % Compute the condition number of matrix Kred-MLM
    cNCLamdaC = cond(Kred-CLambdaC);
    fprintf('>> Matrix Kred-MLM has condition number equal to %d \n\n',cNCLamdaC);
else
    rDLambda = 'undefined';
    cNLambda = 'undefined';
    cNCLamdaC = 'undefined';
    rDCLamdaC = 'undefined';
end

%% 11. Loop over all the fixed point iterations
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Looping over all the fixed point iterations \n');
    fprintf('\t ------------------------------------------- \n');
end
while errordLambda>fixedPointTolerance && counterFixedPointIterations<maxNoIterations
    %% 11i. Update the force vector
    FUpdated = Fred-C*lambdaHat;

    %% 11ii. Solve the structural problem
    dHatred = (Kred-CLambdaC)\FUpdated;
    
    %% 11iii. Solve the perturbed interface coupling equation
    dlambda = -(Lambda\C')*dHatred;
    
    %% 11iv. Compute the updated Lagrange Multiplier field
    lambdaHat = lambdaHat + relaxationFactor*dlambda;
    
    %% 11v. Compute the error in the 2-norm for the increment in the Lagrange Multipliers field
    errordLambda = norm(dlambda);
    
    %% 11vi. Print message on the residual error
    if strcmp(outMsg,'outputEnabled')
        msgFPI = sprintf('\t >> Residual is ||dlambdaHat|| = %d at iteration No. %d \n',errordLambda,counterFixedPointIterations);
        fprintf([reverseStrFixedPointIterations, msgFPI]);
        reverseStrFixedPointIterations = repmat(sprintf('\b'), 1, length(msgFPI));
    end
    
    %% 11vii. Update iteration counter
    counterFixedPointIterations = counterFixedPointIterations+1;
end

%% 12. Print message on the convergence of the fixed point iterations
if strcmp(outMsg,'outputEnabled')
    if counterFixedPointIterations>=maxNoIterations
        fprintf('\t >> Maximum number of iterations has been exceeded \n');
        fprintf('\n');
    else
        fprintf('\t >> Fixed Point iterations converged! \n');
        fprintf('\n');
    end
end

%% 13. Re-assign the output arrays
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

%% 14. Compute the complete force vectors for both patches
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the complete force vectors\n\n');
end
F1 = K1*dHat1;
F2 = K2*dHat2;

%% 15. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Static linear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('____________________________Static Linear Analysis Ended____________________________\n');
    fprintf('####################################################################################\n\n\n');
end

end

