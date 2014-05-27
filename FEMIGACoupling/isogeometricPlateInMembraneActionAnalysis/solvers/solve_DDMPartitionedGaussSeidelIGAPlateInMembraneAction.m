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
function [dHat1,dHat2,lambdaHat,F1,F2] = solve_DDMPartitionedGaussSeidelIGAPlateInMembraneAction(patch1,patch2,...
    lagrangeMultipliers,gaussSeidelTolerance,relaxationFactor,maxNoIterations,intC,outMsg)
%% Function documentation
%
% Returns the displacement field for both patches as well as the discrete
% Lagrange Multiplier field on the Control Points by solving the decomposed
% isogeometric plate in membrane action problem using the partitioned Gauss
% Seidel approach using two methods. Source References :
%
% lagrangeMultipliers.isEnabled == true
% _____________________________________
%
% Ralf Unger et al. "Application of Lagrange multipliers for coupled 
% problems in fluid and structural interactions", Computers & Structures,
% 2007
%
% lagrangeMultipliers.isEnabled == false
% ______________________________________
%
% Tianyang Wang et al. "Concept and Realization of Coupling Software EMPIRE 
% in Multi-Physics Co-Simulation", Computational Methods in Marine 
% Engineering, 2013
%
% It is assumed by default that patch 1 is the master patch and patch 2 is
% the slave patch. Then, patch 1 is always solved as a Neumann problem,
% meaning that it receives forces over the interface and returns a
% displacement field, whereas patch 2 is always solved as a Dirichlet 
% problem, meaning that it receives a displacement field over the interface
% (using a mortar mapping in this algorithm) and returns a new displacement
% as well as force field.
%
%                            Input :
%                    patch1,patch2 : Structures containing all geometrical 
%                                    information (NURBS parameters), 
%                                    technical information (parameters) and 
%                                    the coupling information (structure 
%                                    cb)
%              lagrangeMultipliers : Structure on the Lagrange Multipliers
%                                    discretization as well as on whether
%                                    the interface Lagrange multipliers
%                                    contribute in the convergence of the
%                                    coupled system :
%                                            .field : The Lagrange
%                                                     Multiplier
%                                                     discretization
%                                        .isEnabled : Whether the Lagrange
%                                                     Multiplier field is
%                                                     contributing into the
%                                                     convergence of the
%                                                     solution or not
% isLagrangeMultiplierFieldEnabled : Boolean on whether the Lagrange
%                                    multiplier field is used in the Gauss
%                                    Seidel iterations instead of the
%                                    computation of the coupling matrices
%                                    where it is by default used
%             gaussSeidelTolerance : Convergence tolerance for the fixed 
%                                    point iterations
%                 relaxationFactor : Relaxation factor for the fixed point
%                                    iterations, relaxationFactor \in (0,1]
%                  maxNoIterations : Maximum number of fixed point 
%                                    iterations
%                             intC : Structure containing information on 
%                                    the numerical integration over the 
%                                    coupling surface 
%                           outMsg : Whether or not to output message on 
%                                    refinement progress 
%                                    'outputEnabled' : enables output 
%                                    information
%   
%                           Output :
%                      dHat1,dHat2 : The complete displacement vectors for 
%                                    the two structures
%                        lambdaHat : The Lagrange multipliers field on the 
%                                    interface in a discrete form
%                            F1,F2 : The complete force vector for both 
%                                    patches
%
% Function Layout :
%
% 0. Read input
%
% 1. Compute the stiffness matrices for the two membranes
%
% 2. Determine the direction of the coupling interface
%
% 3. On the choice of the interface Lagrange multipliers field
%
% 4. Compute the coupling matrices for the two patches restricted on the coupling interface
%
% 5. Initialize matrices which are used or updated throughout the Gauss Seidel iterations
%
% 6. Loop over all the Gauss Seidel iterations
% ->
%    6i. (Solve the Neumann Problem) equilibrium equation at patch 1
%
%   6ii. (State Transfer) Transfer the interface displacement field from patch 1 to patch 2
%
%  6iii. (Solve the Dirichlet Problem) equilibrium equation at patch 2 given the displacement field on the interface from patch 1
%
%   6iv. (Stress transfer) Solve for the Lagrange multipliers field on the coupling interface
%
%    6v. Compute the residual of the Lagrange multipliers/displacement field of patch 2 in the 2-norm and check convergence
%
%   6vi. Message on the convergence
%
%  6vii. Update iteration counter
% <-
% 7. Print message on the convergence of the fixed point iterations
%
% 8. Compute the complete force vectors for both patches
%
% 9. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________________________________\n');
    fprintf('###################################################################################\n');
    fprintf('Static linear analysis for the decomposed into two subdomains isogeometric plate in \n');
    fprintf('plane stress action using the partitioned Gauss Seidel method has been initiated \n\n');
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
    fprintf('Gauss Seidel iteration parameters :\n');
    fprintf('-----------------------------------\n');
    fprintf('Iterative method : ');
    if lagrangeMultipliers.isEnabled == false
        fprintf('Tianyang Wang et al. \n');
    elseif lagrangeMultipliers.isEnabled == 1
        fprintf('Ralf Unger et al. \n');
    else
        fprintf(' No iterative method was chosen \n');
        error(' No iterative method was chosen \n');
    end
    fprintf('Gauss Seidel iteration tolerance = %d\n',gaussSeidelTolerance);
    fprintf('Under-relaxation factor = %d\n',relaxationFactor);
    fprintf('Maximum number of Gauss Seidel iterations = %d\n',maxNoIterations);
    fprintf('_________________________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Assign the master and slave relation
isMaster = 1;
isSlave = 0;

% Re-assign all the individual arrays

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
cb1 = patch1.cb;
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
cb2 = patch2.cb;
int2 = patch2.int;

% Number of degrees of freedom
nxi2 = length(CP2(:,1,1));
neta2 = length(CP2(1,:,1));

% Initialize the iteration counter for the Gauss Seidel iterations
counterGaussSeidelIterations = 1;

% Initialize the error
errorMeasure = 2*gaussSeidelTolerance;

% Dimension of the complete displacement vector for patch 1
dHat1 = zeros(length(Fl1),1);

% Get a sequencial numbering of the unconstained DOFs into a vector

% 1st patch :
% ___________

unconstrainedDOFs1 = zeros(length(Fl1),1);
for i=1:length(Fl1)
    unconstrainedDOFs1(i,1) = i;
end
unconstrainedDOFs1(ismember(unconstrainedDOFs1,rb1)) = [];

% 2nd patch :
% ___________

unconstrainedDOFs2 = zeros(length(Fl2),1);
for i=1:length(Fl2)
    unconstrainedDOFs2(i,1) = i;
end
unconstrainedDOFs2(ismember(unconstrainedDOFs2,rb2)) = [];

% Initialize string
if strcmp(outMsg,'outputEnabled')
    reverseStrGSI = '';
end

% On the Lagrange Multiplier field
lagrangeMultiplierField = lagrangeMultipliers.field;

%% 1. Compute the stiffness matrices for the two membranes
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the stiffness matrices for both patches \n');
end

% 1st patch :
% ___________

K1 = computeStiffMtxIGAPlateInMembraneActionLinear(p1,Xi1,q1,Eta1,CP1,isNURBS1,parameters1,int1,'');

% 2nd patch :
% ___________

K2 = computeStiffMtxIGAPlateInMembraneActionLinear(p2,Xi2,q2,Eta2,CP2,isNURBS2,parameters2,int2,'');

%% 2. Determine the direction of the coupling interface
haveTheSameOrientation = findSubdomainInterfaceOrientation(p1,Xi1,q1,Eta1,CP1,isNURBS1,xicoup1,etacoup1,p2,Xi2,q2,Eta2,CP2,isNURBS2,xicoup2,etacoup2);

%% 3. On the choice of the interface Lagrange multipliers field, if selected

if isscalar(lagrangeMultiplierField) && lagrangeMultipliers.isEnabled
   if strcmp(outMsg,'outputEnabled')
        fprintf('>> Automatically choosing the Lagrange multiplier field discretization\n');
    end
elseif ~isscalar(lagrangeMultiplierField) && lagrangeMultipliers.isEnabled
    fprintf('>> Manual choice of the Lagrange multiplier field discretization\n');
end
if isscalar(lagrangeMultiplierField) || lagrangeMultipliers.isEnabled
    % Set the Lagrange mutlipliers variable to NULL
    clear lagrange_multipliers;

    if xicoup2(1)==etacoup2(2)
        % Coupled region in xi-direction for the slave surface
        lagrangeMultiplierField.p = p2;
        lagrangeMultiplierField.Xi = Xi2;

        % Check if we are at the begging or the end of the knot span 
        if etacoup2(1)==Xi2(1)
            lagrangeMultiplierField.CP = zeros(nxi2,length(CP2(1,1,:)));
            lagrangeMultiplierField.CP(:,4) = CP2(:,1,4); 
        else
            % Compute the number of DoFs for the slave patch in the xi-
            % direction
            neta2 = length(CP2(:,1,1));
            lagrangeMultiplierField.CP = zeros(neta2,length(CP2(1,1,:)));
            lagrangeMultiplierField.CP(:,4) = CP2(:,neta2,4); 
        end
    else
        % Coupled region in eta-direction for the slave surface
        lagrangeMultiplierField.p = q2;
        lagrangeMultiplierField.Xi = Eta2;

        % Check if we are at the begging or the end of the knot span 
        if xicoup2(1)==Xi2(1)
            lagrangeMultiplierField.CP = zeros(neta2,length(CP2(1,1,:)));
            lagrangeMultiplierField.CP(:,4) = CP2(1,:,4); 
        else
            % Compute the number of DoFs for the slave patch in the xi-
            % direction
            nxi2 = length(patch2.CP(1,:,1));
            lagrangeMultiplierField.CP = zeros(nxi2,length(CP2(1,1,:)));
            lagrangeMultiplierField.CP(:,4) = CP2(nxi2,:,4); 
        end
    end
end

% Number of DoFs for the Lagrange multipliers
nNodesL = length(lagrangeMultiplierField.Xi)-lagrangeMultiplierField.p-1;
nDOFsL = 2*nNodesL;

% Initial guess for the Lagrange multipliers control variables
lambdaHat = zeros(nDOFsL,1);

%% 4. Compute the coupling matrices for the two patches restricted on the coupling interface
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the Lagrange Multiplier coupling matrices for both patches\n \n');
end

% 1st patch :
% ___________

C1 = computeLagrangeMultiplierCouplingMatrixIGAPlateInMembraneAction(patch1,lagrangeMultiplierField,isMaster,haveTheSameOrientation,intC);

% Compute coupling matrix from equation C1'*d1 + C2'*d2 = 0
CInterface1 = C1(cb1,:)';

% 2nd patch :
% ___________

C2 = computeLagrangeMultiplierCouplingMatrixIGAPlateInMembraneAction(patch2,lagrangeMultiplierField,isSlave,haveTheSameOrientation,intC);

% Compute coupling matrix from equation C1'*d1 + C2'*d2 = 0
CInterface2 = C2(cb2,:)';

% Save C2 as a sparse matrix to accelerate the computations
C2 = sparse(C2);

%% 5. Initialize matrices which are used or updated throughout the Gauss Seidel iterations
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Reducing the stiffness matrix of patch 1 with respect to the constraints\n \n');
end
Kred1 = K1;
Kred1(:,rb1) = [];
Kred1(rb1,:) = [];
Kred1 = sparse(Kred1);

if ~lagrangeMultipliers.isEnabled
    % Initialize the interface force vector applied to patch 1 from patch 2
    FInterface1 = zeros(length(cb1),1);
    
    % Initialize the saved Control Point displacement field for patch 2
    dHatInterface2 = zeros(length(cb2),1);
end

%% 6. Loop over all the Gauss Seidel iterations
if strcmp(outMsg,'outputEnabled')
    fprintf('\t Looping over all the Gauss Seidel iterations \n');
    fprintf('\t -------------------------------------------- \n');
end
while errorMeasure>=gaussSeidelTolerance && counterGaussSeidelIterations<maxNoIterations
    %% 6i. (Solve the Neumann Problem) equilibrium equation at patch 1
    
    % Compute the generalized load vector for the master side
    if lagrangeMultipliers.isEnabled
        % When the Lagrange Multipliers are activated, there is also
        % interface force acting on patch 1 over the interface due to the
        % presence of the multipliers
        FGeneralized1 = Fl1 - C1*lambdaHat;
    else
        % Initialize the force vector acting on patch 1 due to the Neumann
        % conditions as well as the forces over the interface
        FGeneralized1 = Fl1;
        
        % Add the forces resulting from the conservative mapping
        FGeneralized1(cb1) = FGeneralized1(cb1) + FInterface1;
    end
    if norm(FGeneralized1)==0 && counterGaussSeidelIterations==1
        FGeneralized1(length(FGeneralized1),1) = 1;
    end
    
    % Reducing the generalized load vector of patch 1 with respect to the
    % constraints
	Fred1 = FGeneralized1;
    Fred1(rb1) = [];
    
    % Compute the reduced displacement vector
    dHatRed1 = Kred1\Fred1;
    
    % Assemble into the complete displacement vector for patch 1
    dHat1(unconstrainedDOFs1) = dHatRed1;
    
    % Get the interface displacement field for patch 1
    dHatInterface1 = dHat1(cb1);
    
    %% 6ii. (State Transfer) Transfer the interface displacement field from patch 1 to patch 2
    
    % Compute the tranformation matrix using the mortar method
    TMortar = -(CInterface2\CInterface1);
    
    % Compute the interface displacement field on patch 2 via mortar-based
    % mapping (consistent mapping)
    dHatInterfaceTilde2 = TMortar*dHatInterface1;
    
    % Apply a relaxation on the interface displacement field for patch 2 if
    % the Lagrange Multipliers field is not enabled
    if ~lagrangeMultipliers.isEnabled
        % Compute the increment from the previous iteration step
        ddHatInterface2 = dHatInterfaceTilde2 - dHatInterface2;
        
        % Apply a relaxation to the displacement field
        dHatInterface2 = dHatInterface2 + relaxationFactor*ddHatInterface2;
    else
        dHatInterface2 = dHatInterfaceTilde2;
    end
    
    %% 6iii. (Solve the Dirichlet Problem) equilibrium equation at patch 2 given the displacement field on the interface from patch 1
    dHat2 = solveDirichletProblemIGAPlateInMembraneAction(K2,Fl2,rb2,unconstrainedDOFs2,dHatInterface2,cb2);
    
    %% 6iv. (Stress transfer) Solve for the Lagrange multipliers field on the coupling interface
    
    if lagrangeMultipliers.isEnabled
        % Compute the generalized force vector at the current iteration for
        % patch 2
        FGeneralized2Updated = Fl2 - K2*dHat2;
        
        % Something that needs to be done
        FGeneralized2Updated = FGeneralized2Updated(cb2);
        C2Updated = C2(cb2,:);

        % Compute the updated value for the Lagrange multiplier field
        lambdaHatTilde = C2Updated\FGeneralized2Updated;

        % Compute the residual interface forces in terms of the Lagrange
        % multipliers
        dlambda = lambdaHatTilde - lambdaHat;

        % Update the lagrange multipliers field
        lambdaHat = lambdaHat + relaxationFactor*dlambda;
    else        
        % Compute the complete force vector of patch 2
        F2 = K2*dHat2;
        
        % Restrict the force vector over the interface DOFs
        FInterface2 = F2(cb2);
        
        % Compute the force vector applied onto patch 1 over the interface
        % (conservative mapping)
        FInterface1 = -TMortar'*FInterface2;
    end
    %% 6v. Compute the residual of the Lagrange multipliers/displacement field of patch 2 in the 2-norm and check convergence
    if lagrangeMultipliers.isEnabled
        errorMeasure = norm(dlambda);
    else
        errorMeasure = norm(dHatInterface2-dHatInterfaceTilde2);
    end
    
    %% 6vi. Message on the convergence
    if strcmp(outMsg,'outputEnabled')
        msgGSI = sprintf('\t >> Residual is ||dlambdaHat|| = %d at iteration No. %d \n',errorMeasure,counterGaussSeidelIterations);
        fprintf([reverseStrGSI, msgGSI]);
        reverseStrGSI = repmat(sprintf('\b'), 1, length(msgGSI));
    end
    
    %% 6vii. Update iteration counter
    counterGaussSeidelIterations = counterGaussSeidelIterations + 1;
end

%% 7. Print message on the convergence of the fixed point iterations
if strcmp(outMsg,'outputEnabled')
    if counterGaussSeidelIterations>=maxNoIterations
        fprintf('\t >> Maximum number of iterations has been exceeded \n');
        fprintf('\n');
    else
        fprintf('\t >> Gauss Seidel iterations converged! \n');
        fprintf('\n');
    end
end

%% 8. Compute the complete force vectors for both patches
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the complete force vectors\n\n');
end
F1 = K1*dHat1;
F2 = K2*dHat2;

%% 9. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Static linear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('____________________________Static Linear Analysis Ended____________________________\n');
    fprintf('####################################################################################\n\n\n');
end

end

