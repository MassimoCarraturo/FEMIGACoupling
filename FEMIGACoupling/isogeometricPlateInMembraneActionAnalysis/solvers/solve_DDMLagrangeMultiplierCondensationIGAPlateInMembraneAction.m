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
function [dHat1,dHat2,F1,F2,condG] = solve_DDMLagrangeMultiplierCondensationIGAPlateInMembraneAction(patch1,patch2,lagrangeMultiplierField,intC,outMsg)
%% Function documentation
%
% Static condensation of the Lagrange multipliers over the interface for
% the domain decomposition of a plain stress BVP into two subdomains
% within the isogeometric analysis framework.
%
% Gets the geometrical, technical, and coupling information of two
% plates in membrane action together with the constraint and loading 
% conditions and returns the displacement vector for each patch. The number
% of DoFs for the slave patch, i.e. patch 2 must coincide with the number
% of DoFs that discretize the Lagrange multipliers field. 
%
% The algorithm checks if there is any given Lagrange multipliers 
% discretization and uses this, otherwise it creates automatically a 
% discretization for the Lagrange multipliers according to the slave 
% discretization on the coupling interface 
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
%   
%                  Output :
%             dHat1,dHat2 : The complete displacement vectors for the two
%                           patches
%                   F1,F2 : The complete force vectors for both patches
%                   condG : Condition number of the Schur complement matrix
%
% Function Layout :
%
% 0. Read input
%
% 1. Compute the stiffness matrices for the two membranes sorted as K = [Kdd Kdc Kcd Kcc]
%
% 2. Determine the direction of the coupling interface
%
% 3. On the choice of the interface Lagrange multipliers field
%
% 4. Compute the coupling matrices for the two membranes restricted on the coupling interface
%
% 5. Seperate the load vectors into blocks F = [Fd Fc]
%
% 6. Re-arrangement of the numbering of the prescribed DoFs corresponding to u = [uc ud]'
%
% 7. Reduce the stiffness matrix and load vector blocks related to the domains
%
% 8. Compute the Schur complements for both patches
%
% 9. Compute the Schur complement load vectors for both patches
%
% 10. Compute the generalized stiffness matrix for patch 1 (assumed to be the master patch)
%
% 11. Miscellaneous operations
%
% 12. Compute the generalized load vector for the patch 1 (master patch)
%
% 13. Compute the displacement field on the master side of the interface
%
% 14. Compute the displacement field on the slave side of the interface via mortar-based coupling
%
% 15. Compute the domain displacement fields via the Schur projection
%
% 16. Re-assemble the complete displacement vectors
%
% 17. Compute the complete force vectors
%
% 18. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_________________________________________________________________________________________________\n');
    fprintf('################################################################################################\n');
    fprintf('Static linear analysis for the decomposed into two subdomains isogeometric plate in plane stress\n');
    fprintf('action using the Lagrange Multiplier method and static condensation of the Lagrange Multiplier \n');
    fprintf('field with respect to the second patch has been initiated. Master patch is assumed to be patch 1\n');
    fprintf('by default \n\n');
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
    fprintf('Discretization for the Lagrange Multiplier Field equals that \n');
    fprintf('of patch 2 over the coupling interface \n');
    fprintf('_________________________________________________________________________________________________\n\n');

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

% 1st patch :
% ___________

% Number of degrees of freedom
nxi1 = length(CP1(:,1,1));
neta1 = length(CP1(1,:,1));
nDOFS1 = 2*nxi1*neta1;

% Vector containing the numbering of the DOFs in the domain interior 
% Omega\Gamma_c
cDOFsDomain1 = zeros(nDOFS1,1);
for i=1:nDOFS1
    cDOFsDomain1(i) = i;
end
cDOFsDomain1(cb1) = [];

% 2nd patch :
% ___________

% Number of degrees of freedom
nxi2 = length(CP2(:,1,1));
neta2 = length(CP2(1,:,1));
nDOFS2 = 2*nxi2*neta2;

% Vector containing the numbering of the DOFs in the domain interior 
% Omega\Gamma_c
cDOFsDomain2 = zeros(nDOFS2,1);
for i=1:nDOFS2
    cDOFsDomain2(i) = i;
end
cDOFsDomain2(patch2.cb) = [];

%% 1. Compute the stiffness matrices for the two membranes sorted as K = [Kdd Kdc Kcd Kcc]
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the stiffness matrices for both patches and re-arranging into K = [Kdd Kdc; Kcd Kcc] \n');
end

% 1st patch :
% ___________

K1 = computeStiffMtxIGAPlateInMembraneActionLinear(p1,Xi1,q1,Eta1,CP1,isNURBS1,parameters1,int1,'');

% Decomposition and rearrangement of the stiffness matrix
Kdd1 = K1(cDOFsDomain1,cDOFsDomain1);
Kdc1 = K1(cDOFsDomain1,cb1);
% For the Kcd1 use the Kdc1' due to symmetry
Kcc1 = K1(cb1,cb1);

% 2nd patch :
% ___________

K2 = computeStiffMtxIGAPlateInMembraneActionLinear(p2,Xi2,q2,Eta2,CP2,isNURBS2,parameters2,int2,'');

% Decomposition and rearrangement of the stiffness matrix
Kdd2 = K2(cDOFsDomain2,cDOFsDomain2);
Kdc2 = K2(cDOFsDomain2,cb2);
% For the Kcd2 use the Kdc2' due to symmetry
Kcc2 = K2(cb2,cb2);

%% 2. Determine the direction of the coupling interface
haveTheSameOrientation = findSubdomainInterfaceOrientation(p1,Xi1,q1,Eta1,CP1,isNURBS1,xicoup1,etacoup1,p2,Xi2,q2,Eta2,CP2,isNURBS2,xicoup2,etacoup2);

%% 3. On the choice of the interface Lagrange multipliers field
if isscalar(lagrangeMultiplierField)
   if strcmp(outMsg,'outputEnabled')
        fprintf('>> Automatically choosing the Lagrange multiplier field discretization\n');
    end
else
    fprintf('>> Manual choice of the Lagrange multiplier field discretization\n');
end
if isscalar(lagrangeMultiplierField)
    
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

%% 4. Compute the coupling matrices for the two patches restricted on the coupling interface
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the Lagrange Multiplier coupling matrices for both patches\n');
end

% 1st patch :
% ___________

C1 = computeLagrangeMultiplierCouplingMatrixIGAPlateInMembraneAction(patch1,lagrangeMultiplierField,isMaster,haveTheSameOrientation,intC);

% 2nd patch :
% ___________

C2 = computeLagrangeMultiplierCouplingMatrixIGAPlateInMembraneAction(patch2,lagrangeMultiplierField,isSlave,haveTheSameOrientation,intC);

% Restrict the coupling marices on the coupling interface
C1 = C1(cb1,:);
C2 = C2(cb2,:);

%% 5. Seperate the load vectors into blocks F = [Fd Fc]
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Re-arranging the load vectors into blocks F = [Fd; Fc] \n');
end

% 1st patch :
% ___________

% Domain load vector
Fd1 = Fl1(cDOFsDomain1,1);

% Interface load vector
Fc1 = Fl1(cb1,1);

% 2nd patch :
% ___________

% Domain load vector
Fd2 = Fl2(cDOFsDomain2,1);

% Interface load vector
Fc2 = Fl2(cb2,1);

%% 6. Re-arrangement of the numbering of the prescribed DoFs corresponding to u = [uc ud]'
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Re-arranging the DOF numbering to comply with u = [uc; ud]\n');
end

% 1st patch :
% ___________

% Initialiye array
rbcd1 = rb1;

for i=1:length(rb1)
    for j=1:length(cb1) 
        % If the current numbering of a coupling DoF is higher than the 
        % respective for the Dirichlet DoF then decrease the Dirichlet DoF
        % numbering by 1 otherwise make it equal to the initial numbering
        if cb1(j)<=rb1(i)
            rbcd1(i) = rbcd1(i) - 1;
        else
            rbcd1(i) = rbcd1(i);
        end
    end
end

% 2nd patch :
% ___________

% Initialiye array
rbcd2 = rb2;

for i=1:length(rb2)
    for j=1:length(cb2) 
        % If the current numbering of a coupling DoF is higher than the 
        % respective for the Dirichlet DoF then decrease the Dirichlet DoF
        % numbering by 1 otherwise make it equal to the initial numbering
        if cb2(j)<=rb2(i)
            rbcd2(i) = rbcd2(i) - 1;
        else
            rbcd2(i) = rbcd2(i);
        end
    end
end

%% 7. Reduce the stiffness matrix and load vector blocks related to the domains
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Reduce the domain stiffness and coupling matrices and load vectors with respect to contraints\n');
end

% 1st patch :
% ___________

% The block related to the domain is reduced in rows and columns
Kdd1(:,rbcd1) = [];
Kdd1(rbcd1,:) = [];

% The block Kdc is reduced only in the rows
Kdc1(rbcd1,:) = [];

% As block Kcd is used the respective reduced block Kdc' due to symmetry

% The block Fd is reduced in the rows
Fd1(rbcd1) = [];

% 2nd patch :
% ___________

% The block related to the domain is reduced in rows and columns
Kdd2(:,rbcd2) = [];
Kdd2(rbcd2,:) = [];

% The block Kdc is reduced only in the rows
Kdc2(rbcd2,:) = [];

% As block Kcd is used the respective reduced block Kdc' due to symmetry

% The block Fd is reduced in the rows
Fd2(rbcd2) = [];

%% 8. Compute the Schur complements for both patches
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the Schur complements matrices for both patches\n');
end

% 1st patch :
% ___________

invKdd1TimesKdc1 = Kdd1\Kdc1;
SchurComplement1 = Kcc1 - Kdc1'*invKdd1TimesKdc1;

% 2nd patch :
% ___________

invKdd2TimesKdc2 = Kdd2\Kdc2;
SchurComplement2 = Kcc2 - Kdc2'*invKdd2TimesKdc2;

%% 9. Compute the Schur complement load vectors for both patches
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the Schur complement load vectors for both patches\n');
end

% 1st patch :
% ___________

FSchurComplement1 = Fc1 - Kdc1'*(Kdd1\Fd1);

% 2nd patch :
% ___________

FSchurComplement2 = Fc2 - Kdc2'*(Kdd2\Fd2);

%% 10. Compute the generalized stiffness matrix for patch 1 (assumed to be the master patch)
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the generalized stiffness matrix for patch 1\n');
end
G = SchurComplement1 + C1*(C2\SchurComplement2)*(C2'\C1');

%% 11. Miscellaneous operations
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing rank and conditioning of the reduced system matrix \n');

    % size of G:
    [sizeG,~] = size(G);

    % rank of G:
    rankG = rank(G);

    % Conditioning of G
    condG = cond(G);
    fprintf('>> generalized stiffness matrix for master patch: \n');
    fprintf('>> Rank deficiency equal to %d \n',sizeG-rankG);
    fprintf('>> Condition number equal to %d \n',condG);
else
    condG = 'undefined';
end
   

%% 12. Compute the generalized load vector for the patch 1 (master patch)
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the generalized load vector for patch 1 \n');
end
D = FSchurComplement1 - C1*(C2\FSchurComplement2);

%% 13. Compute the displacement field on the master side of the interface
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Solve for the interface displacement field of the patch 1 \n');
end
uHatc1 = G\D;

%% 14. Compute the displacement field on the slave side of the interface via mortar-based coupling
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Solve for the interface displacement field of the patch 2 (mortar-based coupling) \n');
end
uHatc2 = -(C2'\C1')*uHatc1;

%% 15. Compute the domain displacement fields via the Schur projection
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Solve for the domain displacement fields for both patches using the Schur projection\n');
end

% 1st patch :
% ___________

% load vector
FCompleted1 = Fd1 - Kdc1*uHatc1;

% The domain displacement field for the master patch
uHatd1 = Kdd1\FCompleted1;

% 2nd patch :
% ___________

% load vector
FCompleted2 = Fd2 - Kdc2*uHatc2;

% The domain displacement field for the slave patch
uHatd2 = Kdd2\FCompleted2;

%% 16. Re-assemble the complete displacement vectors
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Re-assembling to the complete displacement vector\n');
end

% 1st patch :
% ___________

% Initialize the complete domain displacement vector
uHatdComplete1 = zeros(length(cDOFsDomain1),1);

% Initialize counters
i = 1;
k = 1;
for l = 1:length(cDOFsDomain1)
    if i<=length(rbcd1)
        % if we are in a Dirichlet boundary condition location add 0
        if l==rbcd1(i)
            uHatdComplete1(l,1)=0;
            % update counter
            
            i = i + 1;
        % if not add the Control Point displacement 
        else
            uHatdComplete1(l,1) = uHatd1(k);
            
            % update counter
            k = k + 1;
        end
    else
        uHatdComplete1(l,1) = uHatd1(k);
        
        % update counter
        k = k + 1;
    end
end

% Re-order the complete displacement vector
dHat1 = zeros(length(Fl1),1);

% Coupling DoFs
dHat1(cb1) = uHatc1;

% Domain DoFs
dHat1(cDOFsDomain1) = uHatdComplete1;

% Slave patch:

% Initialize the complete domain displacement vector
uHatdComplete2 = zeros(length(cDOFsDomain2),1);
  
% Initialize counters
i = 1;
k = 1;
for l = 1:length(cDOFsDomain2)
    if i<=length(rbcd2)
        % if we are in a Dirichlet boundary condition location add 0
        if l==rbcd2(i)
            uHatdComplete2(l,1) = 0;
            
            % update counter
            i = i + 1;
        % if not add the Control Point displacement 
        else
            uHatdComplete2(l,1) = uHatd2(k);
            
            % update counter
            k = k + 1;
        end
    else
        uHatdComplete2(l,1) = uHatd2(k);
        
        % update counter
        k = k + 1;
    end
end

% Re-order the complete displacement vector
dHat2 = zeros(length(Fl2),1);

% Coupling DoFs
dHat2(patch2.cb) = uHatc2;

% Domain DoFs
dHat2(cDOFsDomain2) = uHatdComplete2;

%% 17. Compute the complete force vectors
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Computing the complete force vectors\n\n');
end
F1 = K1*dHat1;
F2 = K2*dHat2;

%% 18. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Static linear analysis took %.2d seconds \n\n',computationalTime);
    fprintf('___________________________________Static Linear Analysis Ended___________________________________\n');
    fprintf('##################################################################################################\n\n\n');
end

end

