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
function [dHat1,dHat2,F1,F2,betaTilde,minEig,rankDefinciency,condK] = solve_DDMNitscheIGAPlateInMembraneAction(patch1,patch2,gamma,intC,outMsg)
%% Function documentation
%
% Returns the displacement field corresponding to the application of the
% Nitsche method for the two-domain decomposed isogeometric plate in
% membrane action problem. Source reference :
%
% Apostolatos et al. "A Nitsche-type formulation and comparison of the 
% most common domain decomposition methods in isogeometric analysis",
% International Journal for Numerical Methods in Engineering (2013)
%
%                  Input :
%          patch1,patch2 : Structures containing all geometrical 
%                          information (NURBS parameters), technical 
%                          information (parameters) and the coupling 
%                          information (structure cb)
%                   gamma : The convex linear combination factor
%                    intC : Structure containing information on the
%                           numerical integration over the coupling surface 
%                  outMsg : Whether or not to output message on refinement 
%                           progress 
%                           'outputEnabled' : enables output information
%
%                  Output :
%             dHat1,dHat2 : The complete displacement vectors for the two
%                           patches
%                   F1,F2 : The complete force vectors for both patches
%               betaTilde : The estimated stabilization parameter for the
%                           Nitsche system
%                  minEig : The minimum eigenvalue for the reduced system
%          rankDeficiency : The rank deficiency of the saddle point system
%                   condK : Condition number of the reduced equation system            
%
% Function Layout :
%
% 0. Read input
%
% 1. Compute the stiffness matrices for the two patches
%
% 2. Determine the direction of the coupling interface
%
% 3. Determination of the stabilization factor
%
% 4. Compute the stabilized Nitsche contributions to the stiffness matrices
%
% 5. Compute the stabilized Nitsche contributions to the coupling matrices
%
% 6. Assemble to the global penalty/coupling-stiffness matrix
%
% 7. Compute the global load vector
%
% 8. Reduce the system according to the constraints
%
% 9. Miscellaneous operations
%
% 10. Solve for the global displacement vector
%
% 11. Re-assign the output arrays
%
% 12. Compute the complete force vectors at each patch
%
% 13. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________________________________\n');
    fprintf('#####################################################################################\n');
    fprintf('Static linear analysis for the decomposed into two subdomains isogeometric plate in  \n');
    fprintf('membrane action using the Nitsche-type method [Apostolatos et al.] has been initiated\n\n');
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
    fprintf('Number of elements = %d \n',(length(patch2.Xi) - 2*(patch2.p + 1)+1)*...
        (length(patch2.Eta) - 2*(patch2.q + 1)+1));
    fprintf('_____________________________________________________________________________________\n\n');

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
[rows1,~] = size(K1);
isMaster = true;

% 2nd patch :
% ___________

K2 = computeStiffMtxIGAPlateInMembraneActionLinear(p2,Xi2,q2,Eta2,CP2,isNURBS2,parameters2,int2,'');
[rows2,~] = size(K2);
isSlave = false;

%% 2. Determine the direction of the coupling interface
haveTheSameOrientation = findSubdomainInterfaceOrientation(p1,Xi1,q1,Eta1,CP1,isNURBS1,xicoup1,etacoup1,p2,Xi2,q2,Eta2,CP2,isNURBS2,xicoup2,etacoup2);

%% 3. Determination of the stabilization factor
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the stabilization parameter by solving A*v = lambda B*v \n');
end
betaTilde = computeNitscheStabilParameterIGAPlateInMembraneAction(patch1,K1,patch2,K2,haveTheSameOrientation,intC,outMsg);
if strcmp(outMsg,'outputEnabled')
    fprintf('>> The estimated stabilization parameter is betaTilde = %d \n',betaTilde);
end

%% 4. Compute the stabilized Nitsche contributions to the stiffness matrices
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the stabilized Nitsche contributions to the stiffness matrices\n');
end

% 1st patch :
% ___________

KnStab1 = computeStabilizedNitscheStiffnessMatrixIGAPlateInMembraneAction(patch1,betaTilde,gamma,isMaster,haveTheSameOrientation,intC);
[Kn1,Kp1] = computeStabilizedNitscheStiffnessMatrixIGA4MichaelDebug(patch1,betaTilde,gamma,isMaster,haveTheSameOrientation,intC);

% 2nd patch :
% ___________

KnStab2 = computeStabilizedNitscheStiffnessMatrixIGAPlateInMembraneAction(patch2,betaTilde,gamma,isSlave,haveTheSameOrientation,intC);
[Kn2,Kp2] = computeStabilizedNitscheStiffnessMatrixIGA4MichaelDebug(patch2,betaTilde,gamma,isSlave,haveTheSameOrientation,intC);

%% 5. Compute the stabilized Nitsche contributions to the coupling matrices
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the stabilized Nitsche contributions to the coupling matrices\n');
end

% 1st patch :
% ___________

CnStab1 = computeStabilizedNitscheCouplingMatrixIGAPlateInMembraneAction(patch1,patch2,betaTilde,gamma,isMaster,haveTheSameOrientation,intC);
[Cn,Cp] = computeStabilizedNitscheCouplingMatrixIGA4MichaelDebug(patch1,patch2,betaTilde,gamma,isMaster,haveTheSameOrientation,intC);

% 2nd patch :
% ___________

% Stabilized coupling matrix contribution CnStab2 does not have to be
% computed explicitly since matrix CnStab1' can be used instead due to the
% symmetry of the formulation

%% 6. Assemble to the global penalty/coupling-stiffness matrix
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Assembling to the global Nitsche coupling matrix\n');
end
K = [K1 + KnStab1 CnStab1
     CnStab1'     K2 + KnStab2];
 
%% 7. Compute the global load vector
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Assembling to the global force vector\n');
end
F = [Fl1 
     Fl2];

%% 8. Reduce the system according to the constraints
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Reducing the system according to the given constraints\n');
end

% Rearrange the support vectors with respect to the sequential numbering of
% the full system
for i=1:length(rb2)
    rb2(i) = rows1 + rb2(i);
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

%% 9. Miscellaneous operations
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing minimum eigenvalue, rank and conditioning of the reduced system matrix \n');
    
    % Computing the minimum eigenvalue to the reduced Nitsche system matrix
    [~,eigVal] = eig(Kred);
    
    % Mask the zero values with infty
    eigVal(~eigVal) = Inf;
    minEig = min(min(eigVal));
    
    % Size of the reduced system matrix
    [sizeKred,~] = size(Kred);

    % Rank of the reduced system matrix
    rankKred = rank(Kred);

    % Conditioning of the reduced matrix
    condK = cond(Kred);

    % Message on the minimum eigenvalue of the reduced system matrix
    fprintf('>> The minimum eigenvalue to the reduced system is minEig = %d\n',minEig);
    
    % Message on the rank deficiency of the reduced system matrix
    if sizeKred==rankKred
        rankDefinciency = sizeKred-rankKred;
        fprintf('>> The reduced system has full rank\n');
    else
        rankDefinciency = sizeKred-rankKred;
        fprintf('>> Warning : The reduced system has rank deficiency equal to %d \n',rankDefinciency);
    end
else
    rankDefinciency = 'undefined';
    condK = 'undefined';
    
    % Computing the minimum eigenvalue to the reduced Nitsche system matrix
    [~,eigVal] = eig(Kred);
    
    % Mask the zero values with infty
    eigVal(~eigVal) = Inf;
    minEig = min(min(eigVal));
end

%% 10. Solve for the global displacement vector
Kred = sparse(Kred);
dHatred = Kred\Fred;

%% 11. Re-assign the output arrays
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Re-assembling to the complete displacement vector\n');
end

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
dHat1 = dHat(1:rows1);
dHat2 = dHat(rows1+1:rows1+rows2);

%% 12. Compute the complete force vectors at each patch
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the complete force vectors\n\n');
end
F1 = K1*dHat1;
F2 = K2*dHat2;

%% 13. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Static linear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('___________________Static Linear Analysis Ended____________________\n');
    fprintf('###################################################################\n\n\n');
end

end

