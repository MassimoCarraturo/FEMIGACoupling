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
function patch = fillUpPatch(p,Xi,q,Eta,CP,isNURBS,parameters,Fl,rb,cb,xicoup,etacoup,int)
%% Function documentation
%
% Takes all nessecary arguments that a shell needs and assigns them to
% the structure.
%
%          Input :
%            p,q : The polynomial degrees of the surface
%         Xi,Eta : The knot vectors in xi,eta-directions
%             CP : Set of Control Points and weights
%        isNURBS : Flag on whether the given patch is a NURBS or a B-Spline
%                  patch
%             fl : The load vector
%     parameters : Technical ang geometrical information for the patch
%             rb : fixed DoFs
%             cb : DoFs on coupling interface
% xicoup,etacoup : coupled region
%            int : On the numerical integration
%
%         Output :
%          patch : Structure containing all above information
%
%% Function main body
patch.p = p;
patch.q = q;
patch.Xi = Xi;
patch.Eta = Eta;
patch.CP = CP;
patch.isNURBS = isNURBS;
patch.parameters = parameters;
patch.Fl = Fl;
patch.rb = rb;
patch.cb = cb;
patch.int = int;
patch.xicoup = xicoup;
patch.etacoup = etacoup;

end

