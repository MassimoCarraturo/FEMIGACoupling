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
function [e,minElArea] = computeRelErrorL2CurvedBeamTipShearIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,parameters,internalRadius,externalRadius,pAmp,dHat,error,outMsg)
%% Function documentation
%
% Returns the relative error in the L2 norm for the selected resultant for
% the benchmark case of an isogeometric plate in membrane action problem
% which models a curved beam under tip shear load.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or
%                  B-Spline
%     parameters : Technical and geometrical parameters of the plate
% internalRadius : The internal radius of the curved beam-like plate
% externalRadius : The external radius of the curved beam-like plate
%           pAmp : The applied pressure load in absolute value
%           dHat : The Control Point displacement vector
%          error : Structure related to the resultant error computation
%         outMsg : Whether or not to output message on refinement progress
%                  'outputEnabled' : enables output information
%
%         Output :
%              e : The relative error in the L2-norm for the selected
%                  resultant
%      minElArea : The minimum element area in the IGA mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Get integration rule for the computation of the relative error
%
% 2. Loop over all the elements (knot spans)
% ->
%    2i. Initialize the element area
%
%   2ii. Loop over all the quadrature points
%   ->
%        2ii.1. Compute the Gauss Point coordinates on the NURBS parameter space
%
%        2ii.2. Find the knot span indices
%
%        2ii.3. Compute the NURBS basis functions and/or their derivatives
%
%        2ii.4. Compute the Cartesian image of the current NURBS parametric location
%
%        2ii.5. Compute the curvilinear coordinates of the current NURBS parametric location
%
%        2ii.6. Compute the numerical and the analytical resultant at the parametric location
%
%        2ii.7. Compute the Jacobian of the transformation from the parameter to the physical space and its determinant
%
%        2ii.8. Compute the determinant of the Jacobian to the transformation from the integration to the paramter space
%
%        2ii.9. Compute the element area on the Gauss Point and add the contribution
%
%        2ii.10. Compute the L2-norm of the difference on the Gauss Point and add the contribution
%
%        2ii.11. Compute the L2-norm of the exact resultant on the Gauss Point and add the contribution
%   <-
%  2iii. Check for the minimum element area in the IGA mesh
% <-
% 3. Compute the relative error in the L2-norm
%   
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________________\n');
    fprintf('##################################################################\n');
    fprintf('Computation of the relative error in the L2-norm for the benchmark\n');
    fprintf('problem of a curved beam-like plate in membrane action subject to\n');
    fprintf('tip shear pressure load has been initiated\n\n');
    fprintf('Relative error in the L2-norm for the ');
    if strcmp(error.resultant,'displacement')
        if strcmp(error.component,'x')
            fprintf('displacement component d_x\n');
        elseif strcmp(error.component,'y')
            fprintf('displacement component d_y\n');
        elseif strcmp(error.component,'xy')
            fprintf('displacement vector d\n');
        end
    elseif strcmp(error.resultant,'strain')
        if strcmp(error.component,'x')
            fprintf('strain component epsilon_xx\n');
        elseif strcmp(error.component,'y')
            fprintf('strain component epsilon_yy\n');
        elseif strcmp(error.component,'xy')
            fprintf('strain component epsilon_xy\n');
        end
    elseif strcmp(error.resultant,'stress')
        if strcmp(error.component,'x')
            fprintf('stress component sigma_xx\n');
        elseif strcmp(error.component,'y')
            fprintf('stress component sigma_yy\n');
        elseif strcmp(error.component,'xy')
            fprintf('stress component sigma_xy\n');
        elseif strcmp(error.component,'tensor')
            fprintf('stress tensor sigma\n');
        end
    end
    fprintf('Number of Gauss Points in xi-direction: %d\n',error.xiNGP);
    fprintf('Number of Gauss Points in eta-direction: %d\n',error.etaNGP);
    fprintf('__________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Constant value needed for the computation of the analytical resultants
N = internalRadius^2-externalRadius^2+(internalRadius^2+externalRadius^2)*log(externalRadius/internalRadius);

% Initialize variables
errorInL2 = 0;
exactInL2 = 0;

% Compute the material matrix
D = parameters.E/(1-parameters.nue^2)*[1 parameters.nue 0; parameters.nue 1 0; 0 0 (1-parameters.nue)/2];

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);

% Number of Control Points in xi-,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check the given input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Initialize minimum element area in the IGA mesh by getting some
% indicators from the Control Point net
maxDistanceInxi = 0;
for i=1:neta
    distance = CP(nxi,i,1:3)-CP(1,i,1:3);
    EuclideanDistance = zeros(3,1);
    EuclideanDistance(1:3,1) = distance(1,1,1:3);
    EuclideanDistance = norm(EuclideanDistance);
    if EuclideanDistance>maxDistanceInxi
        maxDistanceInxi = EuclideanDistance;
    end
end
maxDistanceIneta = 0;
for i=1:nxi
    distance = CP(i,neta,1:3)-CP(i,1,1:3);
    EuclideanDistance = zeros(3,1);
    EuclideanDistance(1:3,1) = distance(1,1,1:3);
    EuclideanDistance = norm(EuclideanDistance);
    if EuclideanDistance>maxDistanceIneta
        maxDistanceIneta = EuclideanDistance;
    end
end
minElArea = maxDistanceIneta*maxDistanceIneta;

% Distribute the solution vector into the elements
dHatEl = zeros(mxi-p-1,meta-q-1,2*(p+1)*(q+1));
for j = (q+1):(meta-q-1)
    for i = (p+1):(mxi-p-1)
        % Initialize counter
        k = 1; 
        for c = j-q-1:j-1 
            for b = i-p:i
                dHatEl(i,j,k)   = dHat(2*(c*nxi+b)-1);
                dHatEl(i,j,k+1) = dHat(2*(c*nxi+b));
                
                % update counter
                k = k + 2;
            end
        end
    end
end


%% 1. Get integration rule for the computation of the relative error

% Integration along xi-coordinate line
[xiGP,xiGW] = getGaussPointsAndWeightsOverUnitDomain(error.xiNGP);

% Integration along eta-coordinate line
[etaGP,etaGW] = getGaussPointsAndWeightsOverUnitDomain(error.etaNGP);

%% 2. Loop over all the elements (knot spans)
for j=q+1:length(Eta)-q-1
    for i=p+1:length(Xi)-p-1
        if Eta(j)-Eta(j+1)~=0 && Xi(i)-Xi(i+1)~=0
            %% 2i. Initialize the element area
            elementArea = 0;
            
            %% 2ii. Loop over all the quadrature points
            for l=1:error.etaNGP
                for k=1:error.xiNGP
                    %% 2ii.1. Compute the Gauss Point coordinates on the NURBS parameter space
                    xi = ((1-xiGP(k))*Xi(i)+(1+xiGP(k))*Xi(i+1))/2;
                    eta = ((1-etaGP(l))*Eta(j)+(1+etaGP(l))*Eta(j+1))/2;
                    
                    %% 2ii.2. Find the knot span indices
                    xiSpan = findKnotSpan(xi,Xi,nxi);
                    etaSpan = findKnotSpan(eta,Eta,neta);
                    
                    %% 2ii.3. Compute the NURBS basis functions and/or their derivatives
                    if strcmp(error.resultant,'displacement')
                        nDrv = 0;
                    else
                        nDrv = 1;
                    end
                    dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
                    
                    %% 2ii.4. Compute the Cartesian image of the current NURBS parametric location
                    P = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
                    
                    %% 2ii.5. Compute the curvilinear coordinates of the current NURBS parametric location
                    r = sqrt(P(1,1)^2 + P(2,1)^2);
                    theta = atan(P(2,1)/P(1,1));
                    
                    %% 2ii.6. Compute the numerical and the analytical resultant at the parametric location
                    
                    % Find the Control Point displacement vector affecting
                    % the current knot span
                    dHatActual(:,1) = dHatEl(xiSpan,etaSpan,:);
                    
                    if strcmp(error.resultant,'displacement')
                        error('The displacement error computation for this benchmark has not been implemented');
                    elseif strcmp(error.resultant,'strain')
                        error('The strain error computation for this benchmark has not been implemented');
                    elseif strcmp(error.resultant,'stress')
                        % Compute the Voigt arranged strain tensor
                        epsilonVoigt = computePostprocVoigtStrainIGAPlateInMembraneAction(xiSpan,p,etaSpan,q,CP,dR(:,2:3),dHatActual);
                        
                        % Compute the Voigt arranged stress tensor
                        sigmaVoigt = D*epsilonVoigt;
                        
                        % Compute the analytical stress resultants in the
                        % curvilinear space
                        pr = pAmp*(r+internalRadius^2*externalRadius^2/r^3-(internalRadius^2+externalRadius^2)/r)*sin(theta)/N;
                        pt = pAmp*(3*r-internalRadius^2*externalRadius^2/r^3-(internalRadius^2+externalRadius^2)/r)*sin(theta)/N;
                        prt = -pAmp*(r+internalRadius^2*externalRadius^2/r^3-(internalRadius^2+externalRadius^2)/r)*cos(theta)/N;
                        
                        % Compute the analytical and the numerical stress
                        % resultants over the Cartesian space and at the
                        % Gauss Point location
                        if strcmp(error.component,'x')
                            % Compute the analytical resultant on the Gauss
                            % Point and on the Cartesian space
                            resExact = pr*(cos(theta))^2 + pt*(sin(theta))^2 - 2*prt*sin(theta)*cos(theta);
                            
                            % Compute the numerical resultant on the Gauss
                            % Point 
                            resNumerical = sigmaVoigt(1,1);
                        elseif strcmp(error.component,'y')
                            % Compute the analytical resultant on the Gauss
                            % Point and on the Cartesian space
                            resExact = pr*(sin(theta))^2 + pt*(cos(theta))^2 + 2*prt*sin(theta)*cos(theta);
                            
                            % Compute the numerical resultant on the Gauss
                            % Point
                            resNumerical = sigmaVoigt(2,1);
                        elseif strcmp(error.component,'xy')
                            % Compute the analytical resultant on the Gauss
                            % Point and on the Cartesian space
                            resExact = pr*sin(theta)*cos(theta) - pt*sin(theta)*cos(theta) - prt*((sin(theta))^2-(cos(theta))^2);
                            
                            % Compute the numerical resultant on the Gauss
                            % Point
                            resNumerical = sigmaVoigt(3,1);
                        elseif strcmp(error.component,'tensor')
                            % Compute the analytical resultants on the 
                            % Gauss Point and on the Cartesian space
                            resExact = zeros(2,2);
                            resExact(1,1) = pr*(cos(theta))^2 + pt*(sin(theta))^2 - 2*prt*sin(theta)*cos(theta);
                            resExact(1,2) = pr*sin(theta)*cos(theta) - pt*sin(theta)*cos(theta) - prt*((sin(theta))^2-(cos(theta))^2);
                            resExact(2,1) = resExact(1,2); 
                            resExact(2,2) = pr*(sin(theta))^2 + pt*(cos(theta))^2 + 2*prt*sin(theta)*cos(theta);
                            
                            % Compute the numerical resultants on the Gauss
                            % Point
                            resNumerical = zeros(2,2);
                            resNumerical(1,1) = sigmaVoigt(1,1);
                            resNumerical(1,2) = sigmaVoigt(3,1);
                            resNumerical(2,1) = resNumerical(1,2);
                            resNumerical(2,2) = sigmaVoigt(2,1);
                        end
                    end
                    
                    %% 2ii.7. Compute the Jacobian of the transformation from the parameter to the physical space and its determinant
                    
                    % Initialize Jacobian
                    Jxxi = zeros(2,2);
                    
                    % initialize counter
                    counterJxxi = 0;
                    
                    % Loop over all the non-zero contributions at the span under study
                    for c = 0:q
                        for b = 0:p
                            % Update counter
                            counterJxxi = counterJxxi + 1;
                            % Compute recursively the entries of the Jacobian
                            Jxxi(1,1) = Jxxi(1,1) + CP(xiSpan-p+b,etaSpan-q+c,1)*dR(counterJxxi,2);
                            Jxxi(1,2) = Jxxi(1,2) + CP(xiSpan-p+b,etaSpan-q+c,2)*dR(counterJxxi,2);
                            Jxxi(2,1) = Jxxi(2,1) + CP(xiSpan-p+b,etaSpan-q+c,1)*dR(counterJxxi,3);
                            Jxxi(2,2) = Jxxi(2,2) + CP(xiSpan-p+b,etaSpan-q+c,2)*dR(counterJxxi,3);
                        end
                    end
                    
                    %% 2ii.8. Compute the determinant of the Jacobian to the transformation from the integration to the paramter space
                    detJzetaxi = (Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j))/4;
                    
                    %% 2ii.9. Compute the element area on the Gauss Point and add the contribution
                    elementAreaOnGP = detJzetaxi*det(Jxxi)*xiGW(k)*etaGW(l);
                    elementArea = elementArea + elementAreaOnGP;
                    
                    %% 2ii.10. Compute the L2-norm of the difference on the Gauss Point and add the contribution
                    errorInL2OnGP = norm(resExact-resNumerical)^2*elementAreaOnGP;
                    errorInL2 = errorInL2 + errorInL2OnGP;
                    
                    %% 2ii.11. Compute the L2-norm of the exact resultant on the Gauss Point and add the contribution
                    exactInL2OnGP = norm(resExact)^2*elementAreaOnGP;
                    exactInL2 = exactInL2 + exactInL2OnGP;
                end
            end
            %% 2iii. Check for the minimum element area in the IGA mesh
            if elementArea<minElArea
                minElArea = elementArea;
            end
        end
    end
end
if strcmp(outMsg,'outputEnabled')
    fprintf('>> Minimum element area in the IGA mesh has been found to be %.2d\n',minElArea);
end

%% 3. Compute the relative error in the L2-norm

% Take the square roots of both the diffence and the absolute exact norm
errorInL2 = sqrt(errorInL2);
exactInL2 = sqrt(exactInL2);

% Divide by the norm of the exact solution in the H1
e = errorInL2/exactInL2;
if strcmp(outMsg,'outputEnabled')
    fprintf('>> The relative error in the L2 norm is %d\n\n',e);
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Error computation took %.2d seconds \n\n',computationalTime);
    fprintf('____________________Error Computation Ended_______________________\n');
    fprintf('##################################################################\n\n\n');
end

end

