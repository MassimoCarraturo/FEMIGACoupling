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
function Fl = computeLoadVctAreaIGAPlateInMembraneAction(FlOld,xib,etab,p,q,Xi,Eta,CP,isNURBS,FAmp,direction,int,outMsg)
%% Function documentation
%
% Returns the consistent nodal forces to a area load fload [N/m]. The
% direction of FAmp can be in x, y.
% 
%     Input :
%     FlOld : Existing force vector
%  xib,etab : load extension (e.g. xib = [0 1], etab = [0 1])
%       p,q : the polynomial degrees of the surface
%    Xi,Eta : the knot vectors of the surface
%        CP : The Control Point coordinates and weights
%   isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
%      FAmp : The amplitude of the constant line load or handle to load 
%             function [N/m] for varying loads
% direction : Direction of the applied force:
%             1 = x, 2 = y
%       int : On the numerical integration
%    outMsg : Whether or not to output message on refinement progress
%             'outputEnabled' : enables output information
%
%    Output :
%        Fl : Updated force vector
%
% Function layout :
%
% 0. Read input
%
% 1. On the numerical integration
%
% 2. Find the parametrization of the line integral
%
% 3. Loop over all elements of the load application
% ->
%    3i. Loop over all Gauss points
%    ->
%        3i.1. Compute the coordinate, the map and the Gauss Point location and weight for the fixed parametric coordinate
%
%        3i.2. Compute the IGA basis functions and their first derivatives
%
%        3i.3. Compute the the Jacobian of the transformation from the surface physical to the parameter space and the physical location of the load application
%
%        3i.4. Compute the load amplitude on the Gauss Point is the load is not constant
%
%        3i.5. Compute the product R(xi,eta)*ds in the parameter space or R(xi,eta)*dX in the physical space at the Gauss Point
%
%        3i.6. Compute the determinant of the Jacobian from the physical to the parameter space and the Gauss Point weight
%
%        3i.7. Compute the element load vector
%
%        3i.8. Assemble the element load vector to the global load vector
%    <-
% <-
% 4. Re-assemble the load vector with respect to the global numbering
%
% 5. Update the existing load vector
%
% 6. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('____________________________________________________\n');
    fprintf('####################################################\n');
    if isvector(FlOld)
        fprintf('Update of the load vector corresponding to area\n');
    else
        fprintf('Computation of the load vector corresponding to area\n');
    end
    fprintf('load for the isogeometric plate in membrane action\n');
    fprintf('problem has been initiated\n\n');
    if isnumeric(FAmp)
        fprintf('Constant area load is assumed with amplitude = %.2d\n',FAmp);
    else
        fprintf('Varying area load is assumed\n');
    end
    fprintf('Load extension application Xi x Eta = [%d,%d] x [%d,%d]\n',xib(1),xib(2),etab(1),etab(2));
    fprintf('____________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Length of the knot vectors
mxi = length(Xi);
meta = length(Eta);

% Number of control points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of DOFs
nDOFs = 2*nxi*neta;

% Initialize auxiliary array
Rds = zeros(p+1,q+1);

if isnumeric(FAmp)==1  
    FAmplitude = FAmp;  
end

% Initialize global force vector
F = zeros(nxi,neta,2);

% Initialize element load vector
Fel = zeros(p+1,q+1,2);

% Initialize output array
Fl = zeros(nDOFs,1);

%% 1. On the numerical integration

% Initialize structure
nGP = zeros(2,1);

% Issue the Gauss points for the selected integration scheme:
if strcmp(int.type,'default')
   % Default scheme is the full gaussian quadrature element-wise (FGI)
   nGP(1) = ceil((p+1)/2);
   nGP(2) = ceil((q+1)/2);
elseif strcmp(int.type,'manual')
    % Manual choice of the gauss points
    nGP(1) = int.xiNGPForLoad;
    nGP(2) = int.xetaNGPForLoad;
end

% Issue the Gauss points for the numerical integration
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(nGP(1));
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(nGP(2));

%% 2. Find the parametrization of the surface integral

% For the xi-direction :
% ______________________

% Find the start/end span in xi-direction where to apply the load
i1 = findKnotSpan(xib(1),Xi,nxi);
i2 = findKnotSpan(xib(2),Xi,nxi);   
if xib(2)~=Xi(mxi)  
    i2=i2-1;  
end

% For the eta-direction :
% _______________________

% Find the start/end span in eta-direction where to apply the load
j1 = findKnotSpan(etab(1),Eta,neta);
j2 = findKnotSpan(etab(2),Eta,neta);   
if etab(2)~=Eta(meta)
    j2=j2-1;  
end

%% 3. Loop over all elements of the load application
for j = j1:j2
    for i = i1:i2
        % Check if we are in a non-zero knot span
        if (Xi(i+1)~=Xi(i) && Eta(j+1)~=Eta(j))
            %% 3i. Loop over all Gauss points
            for kv = 1:nGP(2)
                for ku = 1:nGP(1)
                    %% 3i.1. Compute image of the Gauss Points in the NURBS parameter space and the respective Jacobian mappings
                    
                    % Compute the surface parameters
                    xi = (Xi(i+1)+Xi(i) + xiGP(ku)*(Xi(i+1)-Xi(i)))/2;
                    eta = (Eta(j+1)+Eta(j) + etaGP(kv)*(Eta(j+1)-Eta(j)))/2;
                    
                    % Compute the Jacobian determinants for both directions
                    xiMap = (Xi(i+1)-Xi(i))/2;
                    etaMap = (Eta(j+1)-Eta(j))/2;
                    
                    %% 3i.2. Compute the IGA basis functions and their first derivatives
                    nDrv = 1;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface(i,p,xi,Xi,j,q,eta,Eta,CP,isNURBS,nDrv,'');
                    
                    %% 3i.3. Compute the the Jacobian of the transformation from the surface physical to the parameter space and the physical location of the load application
                    
                    % Initialize the Jacobian matrix
                    Jxxi = zeros(2,2);
                    
                    % Initialize the Cartesian coordinates of the load
                    % application point if the load is not constant but
                    % varying
                    if isnumeric(FAmp) == 0
                        x = 0;
                        y = 0;
                    end
                    
                    % Compute the entries of the Jacobian matrix
                    % iteratively
                    counter = 0;
                    for b=0:q
                        for a=0:p
                            % Update counter
                            counter = counter + 1;
                            
                            % Compute the entries of the Jacobian matrix
                            Jxxi(1,1) = Jxxi(1,1) + CP(i-p+a,j-q+b,1)*dR(counter,2);
                            Jxxi(1,2) = Jxxi(1,2) + CP(i-p+a,j-q+b,2)*dR(counter,2);
                            Jxxi(2,1) = Jxxi(2,1) + CP(i-p+a,j-q+b,1)*dR(counter,3);
                            Jxxi(2,2) = Jxxi(2,2) + CP(i-p+a,j-q+b,2)*dR(counter,3);
                            
                            % Compute the physical location of the load
                            % application point if the load is not constant
                            if isnumeric(FAmp) == 0
                                x = dR(counter,1)*CP(i-p+a,j-q+b,1) + x;
                                y = dR(counter,1)*CP(i-p+a,j-q+b,2) + y;
                            end
                        end
                    end 
                    
                    % Compute the determinant of the Jacobian to the
                    % transformation from the physical to the parameter
                    % space
                    detJxxi = det(Jxxi);
                    
                    %% 3i.4. Compute the load amplitude on the Gauss Point is the load is not constant
                    if isnumeric(FAmp) == 0
                        FAmplitude = FAmp(x,y);
                    end
                    
                    %% 3i.5. Compute the product R(xi,eta)*dA in the parameter space or R(xi,eta)*dOmega in the physical space at the Gauss Point
                    
                    % Initialize counter
                    counter = 0;
                    
                    for b = 0:q
                        for a = 0:p
                            % Update counter
                            counter = counter + 1;
                            
                            % Compute the product of the basis function
                            % with the determinant of the Jacobian
                            Rds(a+1,b+1) = dR(counter,1)*detJxxi;
                        end
                    end
                    
                    %% 3i.6. Compute the determinant of the Jacobian from the physical to the parameter space and the Gauss Point weight

                    % Compute the determinant of the Jacobian from the 
                    % physical to the parameter space
                    map = xiMap*etaMap;
                    
                    % Compute the product of the quadrature weights
                    GW = xiGW(ku)*etaGW(kv);
                    
                    %% 3i.7. Compute the element load vector
                    Fel = FAmplitude*map*GW*Rds;
                    
                    %% 3i.8. Assemble the element load vector to the global load vector
                    F(i-p:i,j-q:j,direction) = Fel(:,:) + F(i-p:i,j-q:j,direction);
                end
            end
        end
    end
end

%% 4. Re-assemble the load vector with respect to the global numbering
counter = 1;
for j = 1:length(F(1,:,1))
    for i = 1:length(F(:,1,1))
        % Assemble the x-coordinates of the load vector
        Fl(counter,1) = F(i,j,1);
        
        % Assemble the y-coordinates of the load vector
        Fl(counter+1,1) = F(i,j,2);
        
        % Update counter
        counter = counter + 2;
    end
end

%% 5. Update the existing load vector
if isvector(FlOld)  
    Fl = Fl + FlOld;   
end

%% 6. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    if isvector(FlOld)
        fprintf('Load vector update took %.2d seconds \n\n',computationalTime);
        fprintf('______________Load Vector Update Ended______________\n');
        fprintf('####################################################\n\n\n');
    else
        fprintf('Load vector computation took %.2d seconds \n\n',computationalTime);
        fprintf('____________Load Vector Computation Ended___________\n');
        fprintf('####################################################\n\n\n');
    end
end

end