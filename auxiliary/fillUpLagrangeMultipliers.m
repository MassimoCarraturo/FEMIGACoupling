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
function lm = fillUpLagrangeMultipliers(p,Xi,CP,isNURBS)
%% Function documentation
%
% Takes all nessecary arguments that a membrane needs and assigns them to
% the structure
%
%   Input :
%       p : polynomial degrees of the plane surface
%      Xi : knot vectors in xi,eta-directions
%      CP : set of control points and weights
% isNURBS : Flag on whether the geometrical basis is a NURBS or a B-Spline
%
%  Output :
%      lm : structure containing all above information
%
%% Function main body

lm.p = p;
lm.Xi = Xi;
lm.CP = CP;
lm.isNURBS = isNURBS;

end

