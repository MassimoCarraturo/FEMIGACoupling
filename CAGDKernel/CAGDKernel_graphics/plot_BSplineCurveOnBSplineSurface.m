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
function plot_BSplineCurveOnBSplineSurface(p,q,U,V,CP,isNURBS,grid,u1,v1,u2,v2)
%% Function documentation
%
% Plots the parametric curve on the NURBS surface
%
%    Input :
%      p,q : Polynomial degrees
%      U,V : Knot vectors in u,v-direction
%       CP : Set of control points and weights
%     grid : Number of sampling points to be used
%    u1,v1 : Starting point parameters of the curve
%    u2,v2 : Ending point parameters of the curve
%
%   Output : graphics
%
% Function layout :
%
% 1. Get the coordinates of the sampling points on the curve
%
% 2. Create the geometry
%           
%% Function main body

%% 1. Get the coordinates of the sampling points on the curve
[Xp,Yp,Zp] = createBSplineCurveOnBSplineSurface(p,q,U,V,CP,isNURBS,grid,u1,v1,u2,v2);

%% 2. Create the geometry
line(Xp,Yp,Zp,'Linewidth',2,'color','black');

% hold on;
 
% element edges
% create_el_edges(p,U,CP)

% control points and polygon
%create_conpolygon_curve(CP) 

% axis equal;
% camlight left; lighting phong;
% xlabel('x','FontSize',14);
% ylabel('y','FontSize',14);
% zlabel('z','FontSize',14);
 
% hold off;

end