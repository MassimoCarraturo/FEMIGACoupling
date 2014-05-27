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
function displacement = computeDisplacementVector2D(i,p,u,U,j,q,v,V,CP,d)
%% Function documentation
%
% Returns strain vector displacement = [ux uy]'
%
%        Input : 
%          p,q : polynomial degrees
%          U,V : knot vectors in u,v-direction
%           CP : Control Point coordinates and weights
%          u,v : Parametric coordinate of the NURBS surface
%            d : element displacement vector
%
%       Output :
% displacement : element strain vector strain vector epsilon = [exx eyy 2exy]'
%
% Function layout :
%
% 1. Compute the shape function matrix
% 2. Compute the displacement vector as a product u = N*d
%
%% Function main body

% Number of basis functions per element
ne = (p+1)*(q+1);  

% Number of degrees of freedom per element
ndof = 2*ne;

% Compute NURBS basis functions
R = computeNurbsBasisFunctions2D(i,p,u,U,j,q,v,V,CP);

% Compute the shape function matrix
N = zeros(2,ndof);

for i=1:ne
    N(1,2*i-1) = R(i);
    N(2,2*i) = R(i);
end

% displacement
displacement = N*d;

end