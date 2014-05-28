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
function [N,dN] = computeBilinearBasisFunctionsAndFirstDerivatives(u,v)
%% Function documentation
%
% Returns the bilinear shape functions corresponding to a quadrilateral
% element in a counterclock-wise fashion:
%
%                      eta
%                      / \
%         N4-(-1,+1)    |    N1-(+1,+1)
%              _________|_________
%              |        |        |
%              |        |        |
%              |        |        |
%        -------------------------------> xi
%              |        |        |
%              |        |        |
%              |        |        |
%              _________|_________
%         N3-(-1,-1)    |    N2-(+1,-1)
%                       |
%
%   Input :
%     u,v : The parameters on the parametric domain
%
%  Output :
%       N : A 4x1 array containing the values of the basis functions at the 
%           given parametric locations
%      dN : A 4x2 array containing the values of the derivatives of the
%           basis functions at the given parametric locations as 
%           dN =[dN/dxi; dN/deta]
%
%% Function main body

% Initialize the basis functions array
N = zeros(1,4);

% Initialize the derivatives of the basis functions array
dN = zeros(2,4);

% Assign the values of the functions at the parametric locations
N(1,1) = (1+u)*(1+v)/4;
N(1,2) = (1+u)*(1-v)/4;
N(1,3) = (1-u)*(1-v)/4;
N(1,4) = (1-u)*(1+v)/4;

% Assign the values of the derivatives of the functions at the parametric
% locations

% dN/dxi
dN(1,1) = (1+v)/4;
dN(1,2) = (1-v)/4;
dN(1,3) = -(1-v)/4;
dN(1,4) = -(1+v)/4;

% dN/deta
dN(2,1) = (1+u)/4;
dN(2,2) = -(1+u)/4;
dN(2,3) = -(1-u)/4;
dN(2,4) = (1-u)/4;

end