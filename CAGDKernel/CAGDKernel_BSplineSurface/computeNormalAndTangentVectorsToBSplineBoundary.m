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
function [n,t] = computeNormalAndTangentVectorsToBSplineBoundary(xi,Xi,eta,Eta,GXi,GEta,G3,isOnXi)
%% Function documentation
%
% Returns the normal and the tangent vectors to the shell boundary, given
% the parametric locations of the boundary location given the base vectors. 
% It should be explicitly specified by the last argument of the function in 
% which parametric line the boundary lies on 
%
%    Input : 
%   xi,eta : The parametric locations on the suface boundary only. Interior
%            locations are not allowed and result into an error in the call 
%            of the function
%   Xi,Eta : The knot vectors in u,v-directions
% gXi,gEta : The covariant base vectors on the B-Spline surface
%       G3 : The surface normal base vector
%   isOnXi : The parametric line where the boundary lies on 
%
%   Output :
%        n : The normal to the surface boundary vector
%        t : The tangent to the surface vector (always oriented together 
%            with the tanget to the boundary line base vector)
%
%% Function main body

if isOnXi
    if eta == Eta(1)
        % Tangent vector
        t = GXi/norm(GXi);
        
        % Normal vector
        n = - cross(t,G3);
    elseif eta == Eta(length(Eta))
        % Tangent vector
        t = - GXi/norm(GXi);
        
        % Normal vector
        n = - cross(t,G3);
    end
else
    if xi == Xi(1)
        % Tangent vector
        t = - GEta/norm(GEta);
        
        % Normal vector
        n = - cross(G3,t);
    elseif xi==Xi(length(Xi))
        % Tangent vector
        t = GEta/norm(GEta);
        
        % Normal vector
        n = - cross(G3,t);
    end
end

end

