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
function [dHat1,dHat2,F1,F2,rankDeficiency,condK] = solve_DDMPenaltyIGAPlateInMembraneAction(patch1,patch2,alpha,intC,outMsg)
%% Function documentation
%
% Returns the displacement field and the complete force vectors
% corresponding to the Penalty decomposition method for its application to 
% the isogeometric plate in membrane action analysis decomposed into two 
% subdomains.
%
%                  Input :
%          patch1,patch2 : Structures containing all geometrical 
%                          information (NURBS parameters), technical 
%                          information (parameters) and the coupling 
%                          information (structure cb)
%                   alpha : The penalty parameter
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
%                   F1,F2 : The complete force vectors for both patches
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
% 3. Compute the penalty contributions to the stiffness matrices
%
% 4. Compute the penalty contributions to the coupling matrices
%
% 5. Assemble to the global penalty/coupling-stiffness matrix
%
% 6. Compute the generalized load vector
%
% 7. Reduce the system according to the constraints
%
% 8. Miscellaneous operations
%
% 9. Solve for the global displacement vector
%
% 10. Re-assign the output arrays
%
% 11. Compute the complete force vectors for both patches
%
% 12. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________\n');
    fprintf('#############################################################\n');
    fprintf('Static linear analysis for the decomposed into two subdomains \n');
    fprintf('isogeometric plate in plane stress action using the Penalty \n');
    fprintf('method has been initiated \n\n');
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
    fprintf('Penalty parameter alpha = %d\n',alpha);
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
[rows1,~] = size(K1);
isMaster = 1;

% 2nd patch :
% ___________

K2 = computeStiffMtxIGAPlateInMembraneActionLinear(p2,Xi2,q2,Eta2,CP2,isNURBS2,parameters2,int2,'');
[rows2,~] = size(K2);
isSlave = 0;

%% 2. Determine the direction of the coupling interface
haveTheSameOrientation = findSubdomainInterfaceOrientation(p1,Xi1,q1,Eta1,CP1,isNURBS1,xicoup1,etacoup1,p2,Xi2,q2,Eta2,CP2,isNURBS2,xicoup2,etacoup2);

%% 3. Compute the penalty contributions to the stiffness matrices
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the penalty contributions to the stiffness matrices for both patches\n');
end

% 1st patch :
% ___________

KPenalty1 = computePenaltyStiffnessMatrixIGAPlateInMembraneAction(patch1,alpha,isMaster,haveTheSameOrientation,intC);

% 2nd patch :
% ___________

KPenalty2 = computePenaltyStiffnessMatrixIGAPlateInMembraneAction(patch2,alpha,isSlave,haveTheSameOrientation,intC);

%% 4. Compute the penalty contributions to the coupling matrices
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the penalty contributions to the coupling matrices for both patches\n');
end

% 1st patch :
% ___________

CPenalty1 = computePenaltyCouplingMatrixIGAPlateInMembraneAction(patch1,patch2,alpha,isMaster,haveTheSameOrientation,intC);

% 2nd patch :
% ___________

% No further computation is needed for the penalty contribution to the
% coupling matrix on the slave side cause the transpose of CPenalty1 can be
% used instead due to symmetry

%% 5. Assemble to the global penalty/coupling-stiffness matrix
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Assembling to the global Penalty matrix of the coupled system \n');
end

K = [K1 + KPenalty1 CPenalty1
     CPenalty1'     K2 + KPenalty2];

%% 6. Compute the generalized load vector
F = [Fl1
     Fl2];

%% 7. Reduce the system according to the constraints
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

%% 8. Miscellaneous operations
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing rank and conditioning of the reduced system matrix \n');
    
    % Size of the reduced system matrix
    [sizeKred,~] = size(Kred);

    % Rank of the reduced system matrix
    rankKred = rank(Kred);

    % Conditioning of the reduced matrix
    condK = cond(Kred);
    fprintf('>> The reduced system has condition number equal to %d \n',condK);
    
    % Message on the rank deficiency of the reduced system matrix
    if sizeKred==rankKred
        rankDeficiency = sizeKred-rankKred;
        fprintf('>> The reduced system has full rank\n');
    else
        rankDeficiency = sizeKred-rankKred;
        fprintf('>> Warning : The reduced system has rank deficiency equal to %d \n',rankDeficiency);
    end
else
    rankDeficiency = 'undefined';
    condK = 'undefined';
end

%% 9. Solve for the global displacement vector
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Solving the linear system of %d equations\n',length(Fred));
end
Kred = sparse(Kred);
dHatred = Kred\Fred;

%% 10. Re-assign the output arrays
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

%% 11. Compute the complete force vectors for both patches
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the complete force vectors\n\n');
end
F1 = K1*dHat1;
F2 = K2*dHat2;

%% 12. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Static linear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('_________________Static Linear Analysis Ended_________________\n');
    fprintf('##############################################################\n\n\n');
end

end

