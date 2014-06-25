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
function [K] = computeModifiedLinearSystem...
    (KFEM,KIGA,IBC,dHat,dHatSaved,dHatDot,dHatDotSaved,DOFNumbering,XiB,BSplinePatch,...
    transientAnalysis,t,parameters,int,uFEM,uSavedFEM,uDotFEM,uDotSavedFEM, ...
    DOFNumberingFEM,strMsh,analysisFEM,strDynamics,...
    parametersFEM,Penalty,bodyForces)
%% Function documentation
%
% Returns the linear master stiffness matrix and the complete right-hand 
% side load vector for an isogeometric plate in membrane plane stress 
% action.
%
%             Input : 
%              dHat : Initial guess for the primary field (just an emtpy 
%                     array must be given to this function)
%         dHatSaved : The discrete solution field of the previous time step
%           dHatDot : Initial guess for the time derivative of the primary 
%                     field 
%                     (dummy input for this function)
%      dHatDotSaved : The time derivative of the discrete solution field of 
%                     the previous time step
%                XiB: fixed direction on the interface boundary
%      DOFNumbering : The global numbering of the DOFs in a 3-dimensional 
%                     array
%      BSplinePatch : The polynomial degrees, the knot vectors and the 
%                     Control Point coordinates/weights of the B-Spline 
%                     patch
% transientAnalysis : Transient analysis parameters: (dummy variable to 
%                     this function)
%                               .method : Time integration method
%                            .alphaBeta : Bossak parameter
%                                .gamma : Bossak parameter
%                               .TStart : Start time of the simulation
%                                 .TEnd : End time of the simulation
%                                   .nT : Number of time steps
%                                   .dt : Time step (numeric or adaptive)
%                 t : The current time instance of the transient problem
%        parameters : The parameters of the plane stress problem
%
%            Output :
%                 K : master stiffness matrix
%               
%
% Function layout :
%
% 0. Read input
%
% 1. Choose an integration rule
%
% 1.1 projections
%
% 2. loop over all elements (knot spans)
% ->
%    2i. Compute the determinant of the Jacobian to the transformation from the NURBS space (xi-eta) to the integration domain [-1,1]x[-1,1] 
%
%   2ii. Create an Element Freedom Table
%
%  2iii. Initialize element area size
%
%   2iv. Loop over all Gauss Points
%   ->
%        2iv.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
%
%        2iv.2. Compute the quadrature weight as a tensor product
%
%        2iv.3. Find the correct spans where xi,eta lie in
%
%        2iv.4. Compute the IGA basis functions and their first derivatives at the quadrature point
%
%        2iv.5. Compute the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%        2iv.6. Compute the element area on the Gauss Point
%
%        2iv.7. Compute the element stiffness matrix at the quadrature point multiplying also via the determinants of the Jacobians to the two transformations and the quadrature weight
%
%        2iv.8. Add the contribution from the Gauss Point and assemble to the global matrices/vectors
%   <-
%
%    2v. Check for the minimum element area size in the mesh
% <-
%
%% Function main body

%% 0. Read input

% Get the patch properties
q = BSplinePatch.q;
p = BSplinePatch.p;
Xi = BSplinePatch.Xi;
Eta = BSplinePatch.Eta;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;

% Compute the number of knots and Control Points in both parametric
% directions

neta = length(CP(1,:,1));
nxi = length(CP(:,1,1));

% Number of global DOFs
noDOFs = 2*neta;

% Local number of Control Points
nCPsLoc = (q+1);

% Number of DOFs affected the element under study
noDOFsLoc = 2*(q+1);

% Material matrix D (for plane stress problems)
D = parameters.E/(1-parameters.nue^2)*[1              parameters.nue 0
                                       parameters.nue 1              0 
                                       0              0             (1-parameters.nue)/2];


KpIGA = zeros(noDOFsLoc,noDOFsLoc);

% FEM inputs

% Number of nodes in the mesh
nNodes = length(IBC.nodes(:,1));

% Number of DOFs in the mesh
nDOFsFEM = 2*nNodes;

% Number of nodes at the element level
nNodesElFEM = 3;

% Number of DOFs at the element level
nDOFsElFEM = 2*nNodesElFEM;

% Initializations
RMatrix = zeros(2,noDOFsLoc);
newtonRaphson.eps = 1e-9;
newtonRaphson.maxIt = 100;

% Compute the material matrix for the given problem
if strcmp(analysisFEM.type,'PLANE_STRESS')
    preFactorFEM = parametersFEM.E/(1-parametersFEM.nue^2);
    CFEM = preFactorFEM*[1             parametersFEM.nue 0
                   parametersFEM.nue 1              0
                   0              0             (1-parametersFEM.nue)/2];
elseif strcmp(analysisFEM.physics,'PLANE_STRAIN')
    preFactorFEM = parametersFEM.E*(1-parametersFEM.nue)/(1+parametersFEM.nue)/(1-2*parametersFEM.nue);
    CFEM = preFactorFEM*[1                                 parametersFEM.nue/(1-parametersFEM.nue) 0
                   parametersFEM.nue/(1-parametersFEM.nue) 1                                 0
                   0                                 0                                 (1-2*parametersFEM.nue)/2/(1-parametersFEM.nue)];
end
    



%% 1. Choose an integration rule

% Select the integration scheme
if strcmp(int.type,'default')
    etaNGP = q + 1;
elseif strcmp(int.type,'manual')
    etaNGP = int.etaNGP;
end

% Issue the Gauss Point coordinates and weights
%[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(xiNGP);
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(etaNGP);

% % FEM
if strcmp(int.type,'default')
    nGP = 1;
elseif strcmp(int.type,'manual')
    nGP = int.nGP;
end
[GP,GW] = getGaussRuleOnCanonicalTriangle(nGP);

%% 1.1 Do the projection to find Inserted Knots
P=zeros(1,3); 
projected=zeros(length(IBC.nodes),2); 
for i=1:length(IBC.nodes)
xp=strMsh.nodes(IBC.nodes(i),1);
 yp=strMsh.nodes(IBC.nodes(i),2);
 zp=strMsh.nodes(IBC.nodes(i),3);
P = [xp;yp;zp];
                    xi0 = XiB;
                    eta0 = 0;
                    [xiP,etaP,Projected,flag,noIterations] = ...
                        computeNearestPointProjectionOnBSplineSurface(P,p,Xi,q,Eta,CP,isNURBS,xi0,eta0,newtonRaphson);
                    projected(i,1)=xiP;
                    projected(i,2)=etaP;
end

 Eta_Coupled = sort([Eta' ; projected(:,2)]);
 Eta_Coupled = Eta_Coupled';
 %Xi_Coupled = sort([Xi' ; projected(:,1)]);
 meta = length(Eta);
 

%% 2. loop over all elements (knot spans)
for j = q+1:meta-q-1
  %  loop over sub knot spans of the element
  l = j;
  while Eta_Coupled(l) < Eta(j)
    
        % check if we are in a non-zero knot span
        if  Eta_Coupled(l+1)~=Eta_Coupled(l)
            %% 2i. Compute the determinant of the Jacobian to the transformation from the NURBS space (eta) to the integration domain [-1,1] 
            %
            %         
            %          eta_j+1 - eta_j 
            %  xi,u =  --------------- 
            %                 2 
            %
            
            detJxiu = (Eta_Coupled(l+1)-Eta_Coupled(l))/2;
            
            %% 2ii. Create an Element Freedom Table
            
            % Initialize element freedom table
            EFT = zeros(1,noDOFsLoc);
            
            % Initialize counter
            k = 1;
         
            % Compute the EFT
            for cpj = j-q:j
                    EFT(k)   = DOFNumbering(1,cpj,1);    % DOFNum....(1,cpj,1) first is 1 because 
                    EFT(k+1) = DOFNumbering(1,cpj,2);    % we fix xi at xi=0 and run over eta
                    
                    % update counter
                    k = k + 2;
            end
            
            %% 2iii. Initialize element area size and knot span size
            elLengthSize = 0;
            
            xiSpan = findKnotSpan(Xi(1),Xi,nxi);
            
            %% 2iv. Loop over all Gauss Points
            for keta = 1:length(etaGP)
               
                    %% 2iv.1. Compute the coordinates of the quadrature points at the NURBS domain via mapping
                    xi = XiB;
                    eta = (Eta(j+1)+Eta(j) + etaGP(keta)*(Eta(j+1)-Eta(j)))/2;
                
                    %% 2iv.2. Compute the quadrature weight as a tensor product (check if this is right! We avoid the tensor product...)
                    GW = etaGW(keta);
                
                    %% 2iv.3. Find the correct spans where,eta lie in
                    etaSpan = findKnotSpan(eta,Eta,neta);
                
                    %% 2iv.4. Compute the IGA basis functions and their first derivatives at the quadrature point
                    nDrv = 1;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
                    
                    %% 2iv.5. Compute the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta) and the physical coordinates of the Gauss point
                
                    % Initialize Jacobian
                    Jxxi = zeros(1,2);
                    
                    % Initialize Cartesian coordinates
                    x = 0;
                    y = 0;
                    z = 0;
                
                    % initialize counter
                    k = 0;
                
                    % Loop over all the non-zero contributions at the span
                    % under study
                    for c = 0:q
                        for b = 0:p
                            % Update counter
                            k = k + 1;
                        
                            % Compute recursively the entries of the Jacobian
%                             Jxxi(1,1) = Jxxi(1,1) + CP(p,j-q+c,1)*dR(k,2);
%                             Jxxi(1,2) = Jxxi(1,2) + CP(p,j-q+c,2)*dR(k,2);
                            Jxxi(1,1) = Jxxi(1,1) + CP(xiSpan-p+b,j-q+c,1)*dR(k,3);
                            Jxxi(1,2) = Jxxi(1,2) + CP(xiSpan-p+b,j-q+c,2)*dR(k,3);
                            
                            % Compute the Cartesian coordinates of the
                            % Gauss Point
                            x = x + CP(xiSpan-p+b,j-q+c,1)*dR(k,1);
                            y = y + CP(xiSpan-p+b,j-q+c,2)*dR(k,1);
                            z = z + CP(xiSpan-p+b,j-q+c,3)*dR(k,1);
                            a=[x y z];
                        end
                    end
              
                
                    % Compute the determinant of the Jacobian
                    detJxxi = norm(Jxxi);
                    
                    
                    %% 2iv.6. Compute the element length on the Gauss Point
                    elLengthSizeOnGP = abs(detJxxi)*abs(detJxiu)*GW;
                    elLengthSize = elLengthSize + elLengthSizeOnGP;
                
                    %% 2iv.7. 
                    for counter = 1:nCPsLoc
                        RMatrix(1,2*counter-1) = dR(counter,1);
                        RMatrix(2,2*counter) = dR(counter,1);
                    end
                    
%                     b = bodyForces(x,y,z,t);
%                     [Ke,Fe] = computeElStiffMtxAndLoadVctPlateInMembraneActionLinear(nCPsLoc,dR(:,1),dR(:,2:3),Jxxi,D,b(1:2,1));
                    
                    %% 2iv.8. Add the contribution from the Gauss Point and assemble to the global matrices/vectors
                    
                    % For the stiffness matrix
                    KpIGA = Penalty*(RMatrix'*RMatrix)*elLengthSizeOnGP;
                    KIGA(EFT, EFT) = KIGA(EFT, EFT) + KpIGA;   
            end
%             %% 2v. Check for the minimum element area size in the mesh
%             if elLengthSize < minElSize
%                 minElSize = elLengthSize;
%             end
        end
    l = l+1;
  end
end

% FEM
for counterEl = 1:length(strMsh.elements(:,1))
    %% 2i. Get the element in the mesh
    element = strMsh.elements(counterEl,:);
    
    %% 2ii. Get the nodes in the element
    node1 = strMsh.nodes(element(1,1),:);
    node2 = strMsh.nodes(element(1,2),:);
    node3 = strMsh.nodes(element(1,3),:);
    
    %% 2iii. Create an Element Freedome Table (EFT)
    EFT = zeros(nDOFsElFEM,1);
    for counterEFT = 1:nNodesElFEM
        EFT(2*counterEFT-1) = 2*element(1,counterEFT)-1;
        EFT(2*counterEFT) = 2*element(1,counterEFT);
    end
    
    %% 2iv. Loop over all the quadrature points
    for counterGP = 1:nGP
        %% 2iv.1. Transform the Gauss Point location from the parameter to the physical space
        xGP = GP(counterGP,1)*node1(1,:) + GP(counterGP,2)*node2(1,:) + (1-GP(counterGP,1)-GP(counterGP,2))*node3(1,:);
        
        %% 2iv.2. Compute the basis functions and their derivatives at the Gauss Point
        [dN,detJxxi] = computeCST2DBasisFunctionsAndFirstDerivatives(node1,node2,node3,xGP(1,1),xGP(1,2));
        
        %% 2iv.3. Form the basis functions matrix at the Gauss Point
        N = [dN(1,1) 0       dN(2,1) 0       dN(3,1) 0
             0       dN(1,1) 0       dN(2,1) 0       dN(3,1)];
        
        %% 2iv.4. Form the B-Operator matrix for the plate in membrane action problem
        B = [dN(1,2) 0       dN(2,2) 0       dN(3,2) 0
             0       dN(1,3) 0       dN(2,3) 0       dN(3,3)
             dN(1,3) dN(1,2) dN(2,3) dN(2,2) dN(3,3) dN(3,2)];
         
        %% 2iv.5. Compute the element stiffness matrix at the Gauss Point and assemble to master stiffness matrix via the EFT
        KpFEM = Penalty*(N'*N)*elLengthSizeOnGP;
        KFEM(EFT, EFT) = KFEM(EFT, EFT) + KpFEM;   
        
        
%         %% 2iv.6. Compute the element load vector due to body forces and assemble to master load vector via the EFT
%         bF = bodyForces(xGP(1,1),xGP(1,2),xGP(1,3));
%         FFEM(EFT) = FFEM(EFT) + N'*bF(1:2,1)*detJxxi*GW(counterGP);
    end
end


%% 5. Generate Sparse Stiffness Matrix
K=blkdiag(KFEM,KIGA);

end