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
function Fl = computeLoadVctPointIGAPlateInMembraneAction(FlOld,xi,eta,p,q,Xi,Eta,CP,isNURBS,FAmp,direction,outMsg)
%% Function documentation
%
% Returns the consistent load vector corresponding to the application of a
% point load at the parametric location for the isogeometric plate in
% membrane action.
%
%       Input :
%       FlOld : Existing force vector
%      xi,eta : Parametric location of the load application
%         p,q : The polynomial degrees in xi-,eta- direction
%      Xi,Eta : The knot vectors in xi-,eta- direction
%          CP : The Control Point coordinates and weights
%     isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
%        FAmp : The amplitude of the point load
%   direction : The direction of the point load:
%               1 = x, 2 = y
%      outMsg : Whether or not to output message on refinement progress
%               'outputEnabled' : enables output information
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the IGA basis functions at the load application parametric location
%
% 2. Compute the consistent load vector corresponding to the application of a point load
%
% 3. Re-assemble the load vector with respect to the global numbering
%
% 4. Update the existing load vector
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('_____________________________________________________\n');
    fprintf('#####################################################\n');
    if isvector(FlOld)
        fprintf('Update of the load vector corresponding to point\n');
    else
        fprintf('Computation of the load vector corresponding to point\n');
    end
    fprintf('boundary load for the isogeometric plate in membrane\n');
    fprintf('action problem has been initiated\n\n');
    fprintf('Point load application at surface parameters (%d,%d)\n',xi,eta);
    fprintf('_____________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of Control Points in xi-, eta- direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of Control Points
nNodes = nxi*neta;

% Number of DOFs
nDOFs = 2*nNodes;

% Find the span where the point load is applied
xiSpan = findKnotSpan(xi,Xi,nxi);
etaSpan = findKnotSpan(eta,Eta,neta);

% Initialize the load vector
F = zeros(nxi,neta,2);

% Initialize output array
Fl = zeros(nDOFs,1);

%% 1. Compute the IGA basis functions at the load application parametric location
nDrv = 0;
R = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv,'');

%% 2. Compute the consistent load vector corresponding to the application of a point load

% initialize counter
k = 1;

% Loop over all the basis function contributions
for c = 0:q
    for b = 0:p
        % Sum up the contributions from aech basis function
        F(xiSpan-p+b,etaSpan-q+c,direction) = R(k,1)*FAmp;
        
        % Update counter
        k = k + 1;
    end
end

%% 3. Re-assemble the load vector with respect to the global numbering
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

%% 4. Update the existing load vector
if isvector(FlOld)  
    Fl = Fl + FlOld;   
end

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    if isvector(FlOld)
        fprintf('Load vector update took %.2d seconds \n\n',computationalTime);
        fprintf('______________Load Vector Update Ended_______________\n');
        fprintf('#####################################################\n\n\n');
    else
        fprintf('Load vector computation took %.2d seconds \n\n',computationalTime);
        fprintf('____________Load Vector Computation Ended____________\n');
        fprintf('#####################################################\n\n\n');
    end
end

end