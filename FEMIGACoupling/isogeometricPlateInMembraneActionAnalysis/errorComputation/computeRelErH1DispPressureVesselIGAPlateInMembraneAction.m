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
function [e,minElArea] = computeRelErH1DispPressureVesselIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,parameters,internalRadius,externalRadius,pAmp,dHat,error,outMsg)
%% Function documentation
%
% Returns the relative error in the H1 norm for the selected displacement
% component for the benchmark case of an isogeometric plate in membrane 
% action problem which models the one quarter of a pressure vessel plate 
% subject to uniform internal pressure.
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
%              e : The relative error in the H1-norm for the selected
%                  displacement component
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
%    ->
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
%        2ii.6. Compute the displacement vector u=[ux uy]' at the parametric location
%
%        2ii.7. Compute the Jacobian of the transformation from the parameter to the physical space and its determinant
%
%        2ii.8. Compute the displacement gradient tensor gradient(u) = [ux,x ux,y ; uy,x uy,y]'
%
%        2ii.9. Compute the numerical and the analytical resultant at the parametric location
%
%        2ii.10. Compute the determinant of the Jacobian to the transformation from the integration to the paramter space
%
%        2ii.11. Compute the element area on the Gauss Point and add the contribution
%
%        2ii.12. Compute the H1-norm of the difference on the Gauss Point and add the contribution
%
%        2ii.13. Compute the H1-norm of the exact resultant on the Gauss Point and add the contribution
%    <-
%  2iii. Check for the minimum element area in the IGA mesh
% <-
% 3. Compute the relative error in the H1-norm
%   
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________________________\n');
    fprintf('#####################################################################\n');
    fprintf('Computation of the relative error in the H1-norm on the displacements\n');
    fprintf('for the benchmark problem of a infinite plate with a circular hole in\n'); 
    fprintf('membrane action subject to tensional loading in x=+-infty has been\n');
    fprintf('initiated\n\n');
    fprintf('Relative error in the H1-norm for the ');
    if strcmp(error.component,'x')
        fprintf('displacement component d_x\n');
    elseif strcmp(error.component,'y')
        fprintf('displacement component d_y\n');
    elseif strcmp(error.component,'tensor')
        fprintf('displacement vector d\n');
    end
    fprintf('Number of Gauss Points in xi-direction: %d\n',error.xiNGP);
    fprintf('Number of Gauss Points in eta-direction: %d\n',error.etaNGP);
    fprintf('_____________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Initialize variables
errorInH1 = 0;
exactInH1 = 0;

% Number of knots in xi-,eta-direction
mxi = length(Xi);
meta = length(Eta);

% Number of Control Points in xi-,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check the given input
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Local number of Control Points
nNodeLoc = (p+1)*(q+1);

% Local number of DOFs
nDOFsLoc = 2*nNodeLoc;

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
                    nDrv = 1;
                    dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv,'');
                    
                    %% 2ii.4. Compute the Cartesian image of the current NURBS parametric location
                    P = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));
                    
                    %% 2ii.5. Compute the curvilinear coordinates of the current NURBS parametric location
                    r = sqrt(P(1,1)^2 + P(2,1)^2);
                    theta = 2*pi() - atan(P(2,1)/(-P(1,1)));
                    
                    %% 2ii.6. Compute the displacement vector u=[ux uy]' at the parametric location
                    
                    % Find the Control Point displacement vector affecting
                    % the current knot span
                    dHatActual(:,1) = dHatEl(xiSpan,etaSpan,:);
                    
                    % Compute the displacement vector u=[ux uy]'
                    d = computePostprocDisplacementIGAPlateInMembraneAction(p,q,dR(:,1),dHatActual);
                    
                    %% 2ii.7. Compute the Jacobian of the transformation from the parameter to the physical space and its determinant
                    
                    % Initialize Jacobian
                    Jxxi = zeros(2,2);
                    
                    % Initialize counter
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
                    
                    %% 2ii.8. Compute the displacement gradient tensor gradient(u) = [ux,x ux,y ; uy,x uy,y]'
                    
                    % Initialize matrix containing the derivatives for each 
                    % basis function at the Cartesian space
                    dRdx = zeros(nNodeLoc,2);
                    
                    % Initialize the B-operator matrices for both 
                    % displacement components
                    uxBGradient = zeros(2,nDOFsLoc);
                    uyBGradient = zeros(2,nDOFsLoc);
                    
                    % Initialize matrix containing the displacement
                    % gradient components
                    displacementGradient = zeros(2,2);
                    
                    for c = 1:nNodeLoc
                        % Compute the derivatives of the basis
                        % functions in the Cartesian space
                        dRdx(c,:) = Jxxi\dR(c,2:3)';
                        
                        % Compute the entries of the B-operator
                        % matrices for the components of the
                        % displacement gradient tensor
                        uxBGradient(1,2*c-1) = dRdx(c,1);
                        uxBGradient(2,2*c-1) = dRdx(c,2);
                        uyBGradient(1,2*c) = dRdx(c,1);
                        uyBGradient(2,2*c) = dRdx(c,2);
                    end
                    
                    % Compute the displacement gradient sorted in a matrix
                    uxGradient = uxBGradient*dHatActual;
                    uyGradient = uyBGradient*dHatActual;
                    
                    % Compute displacement gradient tensor in a matrix form
                    % gradient(u) = [ux,x ux,y; uy,x uy,y]'
                    displacementGradient(1,1) = uxGradient(1);
                    displacementGradient(1,2) = uxGradient(2);
                    displacementGradient(2,1) = uyGradient(1);
                    displacementGradient(2,2) = uyGradient(2);
                    
                    %% 2ii.9. Compute the numerical and the analytical resultant at the parametric location
                    
                    % Compute the exact displacement components over the
                    % curvlinear basis
                    ur = (1+parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)*r/(1+parameters.nue)+externalRadius^2/r)/parameters.E/(externalRadius^2-internalRadius^2);
                    ut = 0;

                    % Compute the analytical and the numerical stress
                    % resultants over the Cartesian space and at the
                    % Gauss Point location
                    if strcmp(error.component,'x')
                        % Compute the analytical resultant on the Gauss
                        % Point and on the Cartesian space:
                        
                        % The L2-part of the norm
                        resExactL2 = ur*cos(theta)-r*ut*sin(theta);

                        % The H1-part of the norm
                        resExactH1 = zeros(2,1);
                        resExactH1(1,1) = (1+parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)/...
                            (1+parameters.nue)-externalRadius^2/r^2)*cos(theta)^2/(parameters.E *(externalRadius^2-...
                            internalRadius^2))+sin(theta)^2*(1+parameters.nue)*internalRadius^2*pAmp*((1-...
                            2*parameters.nue)*r/(1+parameters.nue)+externalRadius^2/r)/...
                            (r*parameters.E *(externalRadius^2-internalRadius^2));
                        resExactH1(2,1) = cos(theta)*(1+parameters.nue)*internalRadius^2*pAmp*((1-...
                            2*parameters.nue)/(1+parameters.nue)-externalRadius^2/r^2)*sin(theta)/...
                            (parameters.E *(externalRadius^2-internalRadius^2))-cos(theta)*(1+...
                            parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)*r/(1+parameters.nue)+...
                            externalRadius^2/r)*sin(theta)/(r*parameters.E *(externalRadius^2-internalRadius^2));
                        
                        % Compute the numerical resultant on the Gauss
                        % Point:
                        
                        % The L2-part of the norm
                        resNumericalL2 = d(1,1);
                        
                        % The H1-part of the norm
                        resNumericalH1 = [displacementGradient(1,1)
                                          displacementGradient(1,2)];
                    elseif strcmp(error.component,'y')
                        % Compute the analytical resultant on the Gauss
                        % Point and on the Cartesian space:
                        
                        % The L2-part of the norm
                        resExactL2 = ur*sin(theta)+r*ut*cos(theta);

                        % The H1-part of the norm
                        resExactH1 = zeros(2,1);
                        resExactH1(1,1) = cos(theta)*(1+parameters.nue)*internalRadius^2*pAmp*((1-...
                            2*parameters.nue)/(1+parameters.nue)-externalRadius^2/r^2)*sin(theta)/...
                            (parameters.E*(externalRadius^2-internalRadius^2))-cos(theta)*(1+...
                            parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)*r/...
                            (1+parameters.nue)+externalRadius^2/r)*sin(theta)/...
                            (r*parameters.E*(externalRadius^2-internalRadius^2));
                        resExactH1(2,1) = (1+parameters.nue)*internalRadius^2*pAmp*((1-...
                            2*parameters.nue)/(1+parameters.nue)-externalRadius^2/r^2)*sin(theta)^2/...
                            (parameters.E*(externalRadius^2-internalRadius^2))+cos(theta)^2*(1+...
                            parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)*r/(1+parameters.nue)+...
                            externalRadius^2/r)/(r*parameters.E*(externalRadius^2-internalRadius^2));
                        
                        % Compute the numerical resultant on the Gauss
                        % Point
                        
                        % The L2-part of the norm
                        resNumericalL2 = d(2,1);
                        
                        % The H1-part of the norm
                        resNumericalH1 = [displacementGradient(2,1)
                                          displacementGradient(2,2)];
                    elseif strcmp(error.component,'tensor')
                        % Compute the analytical resultants on the 
                        % Gauss Point and on the Cartesian space
                        
                        % The L2-part of the norm
                        dxExact = ur*cos(theta)-r*ut*sin(theta);
                        dyExact = ur*sin(theta)+r*ut*cos(theta);
                        resExactL2 = zeros(2,1);
                        resExactL2(1,1) = dxExact;
                        resExactL2(2,1) = dyExact;
                        
                        % The H1-part of the norm
                        resExactH1 = zeros(2,2);
                        resExactH1(1,1) = (1+parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)/...
                            (1+parameters.nue)-externalRadius^2/r^2)*cos(theta)^2/(parameters.E*(externalRadius^2-...
                            internalRadius^2))+sin(theta)^2*(1+parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)*r/...
                            (1+parameters.nue)+externalRadius^2/r)/(r*parameters.E*(externalRadius^2-internalRadius^2));
                        resExactH1(1,2) = cos(theta)*(1+parameters.nue)*internalRadius^2*pAmp*((1-...
                            2*parameters.nue)/(1+parameters.nue)-externalRadius^2/r^2)*sin(theta)/...
                            (parameters.E*(externalRadius^2-internalRadius^2))-cos(theta)*(1+...
                            parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)*r/(1+...
                            parameters.nue)+externalRadius^2/r)*sin(theta)/...
                            (r*parameters.E*(externalRadius^2-internalRadius^2));
                        resExactH1(2,1) = cos(theta)*(1+parameters.nue)*internalRadius^2*pAmp*((1-...
                            2*parameters.nue)/(1+parameters.nue)-externalRadius^2/r^2)*sin(theta)/...
                            (parameters.E*(externalRadius^2-internalRadius^2))-cos(theta)*(1+...
                            parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)*r/(1+...
                            parameters.nue)+externalRadius^2/r)*sin(theta)/...
                            (r*parameters.E*(externalRadius^2-internalRadius^2));
                        resExactH1(2,2) = (1+parameters.nue)*internalRadius^2*pAmp*((1-2*parameters.nue)/...
                            (1+parameters.nue)-externalRadius^2/r^2)*sin(theta)^2/(parameters.E*(externalRadius^2-...
                            internalRadius^2))+cos(theta)^2*(1+parameters.nue)*internalRadius^2*pAmp*((1-...
                            2*parameters.nue)*r/(1+parameters.nue)+externalRadius^2/r)/...
                            (r*parameters.E*(externalRadius^2-internalRadius^2));
                        
                        % Compute the numerical resultants on the Gauss
                        % Point
                        
                        % The L2-part of the norm
                        resNumericalL2 = d;
                        
                        % The H1-part of the norm
                        resNumericalH1 = displacementGradient;
                    end
                    
                    %% 2ii.10. Compute the determinant of the Jacobian to the transformation from the integration to the paramter space
                    detJzetaxi = (Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j))/4;
                    
                    %% 2ii.11. Compute the element area on the Gauss Point and add the contribution
                    elementAreaOnGP = detJzetaxi*det(Jxxi)*xiGW(k)*etaGW(l);
                    elementArea = elementArea + elementAreaOnGP;
                    
                    %% 2ii.12. Compute the H1-norm of the difference on the Gauss Point and add the contribution
                    % If the relative error computation is performed for a
                    % tensor, then it is needed the double contraction of
                    % the difference to the tensors for the computation of
                    % the H1-norm
                    errorInH1OnGP = (norm(resExactL2-resNumericalL2)^2 + norm(resExactH1-resNumericalH1)^2)*elementAreaOnGP;
                    errorInH1 = errorInH1 + errorInH1OnGP;
                    
                    %% 2ii.13. Compute the H1-norm of the exact resultant on the Gauss Point and add the contribution
                    % If the relative error computation is performed for a
                    % tensor, then it is needed the double contraction of
                    % the tensors for the computation of the H1-norm
                    exactInH1OnGP = (norm(resExactL2)^2 + norm(resExactH1)^2)*elementAreaOnGP;
                    exactInH1 = exactInH1 + exactInH1OnGP;
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

%% 3. Compute the relative error in the H1-norm

% Take the square roots of both the diffence and the absolute exact norm
errorInH1 = sqrt(errorInH1);
exactInH1 = sqrt(exactInH1);

% Divide by the norm of the exact solution in the H1
e = errorInH1/exactInH1;
if strcmp(outMsg,'outputEnabled')
    fprintf('>> The relative error in the H1 norm is %d\n\n',e);
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Error computation took %.2d seconds \n\n',computationalTime);
    fprintf('______________________Error Computation Ended________________________\n');
    fprintf('#####################################################################\n\n\n');
end

end

