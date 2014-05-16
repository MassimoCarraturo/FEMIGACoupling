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
function N = computeBilinearBasisFunctions(u,v)
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
%
%% Function main body

% Initialize the basis functions array
N = zeros(1,4);

% Assign the values at the parametric locations
N(1,1) = (1+u)*(1+v)/4;
N(1,2) = (1+u)*(1-v)/4;
N(1,3) = (1-u)*(1-v)/4;
N(1,4) = (1-u)*(1+v)/4;

end