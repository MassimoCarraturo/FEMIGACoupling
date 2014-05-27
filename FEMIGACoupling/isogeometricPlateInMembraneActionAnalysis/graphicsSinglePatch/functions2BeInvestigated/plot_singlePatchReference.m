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
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = plot_singlePatchReference(p,q,U,V,CP,rb,fl,graph)
%% Function documentation
% Plots two NURBS patches in the reference configuration
%
%         Input : 
%           p,q : polynomial degrees of the NURBS surface
%           U,V : knot vectors in U,V-directions
%            CP : set of control points and weights
%            rb : Structure containing information on the supports
%            fl : The force vector
%         graph : structure containing all information on the plots
%
%        Output :
%         index : the index of the graph
%
%% Function main body

% Booleans
IsUndeformed = 0;

% Number the current plot
figure(graph.index)

% Plot the NURBS patch
[Xp,Yp,Zp] = create_surface(p,q,U,V,CP,50,50);
surf(Xp,Yp,Zp,'FaceColor','green','EdgeColor','none');
hold;
create_edges(p,q,U,V,CP,IsUndeformed,50,50)

% control points and control polygon
create_control_polygon(CP)

% Create the supports
[xs,ys,zs] = create_supports(CP,rb);

% Plot the supports
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end

% Create force arrows
[xf,yf,zf] = create_force_arrows(CP,fl);
    
% Plot the load arrows
for k =1:length(xf(:,1))
    plot3(xf(k,:),yf(k,:),zf(k,:),'Linewidth',5);
    plot3(xf(k,1),yf(k,1),zf(k,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
end

% Adjust graph settings
camlight left; 
lighting phong;
view(2);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title ('Reference configuration');
hold off;

% Update plot index by 1
index = graph.index+1;

end

