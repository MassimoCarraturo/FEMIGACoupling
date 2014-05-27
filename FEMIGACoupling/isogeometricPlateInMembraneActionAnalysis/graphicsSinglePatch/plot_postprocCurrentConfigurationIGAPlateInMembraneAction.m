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
function  plot_postprocCurrentConfigurationIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,xiGrid,etaGrid,rb,Fl,dHat,graph)
%% Function documentation
%
% Plots a window with the reference and/or the current configuration of an
% IGA plate in membrane action.
%
%          Input :
%            p,q : Polynomial degrees
%         Xi,Eta : Knot vectors in xi,eta-direction
%             CP : Control point coordinates and weights of the undeformed 
%                  plate
%        isNURBS : Flag on whether the geometrical basis is NURBS or 
%                  B-Spline
% xiGrid,etaGrid : The grid points used for the plotting of the NURBS
%                  geometry
%             rb : Vector containing information on the supports
%             Fl : The applied load vector
%           dHat : The displacement field of the control points
%          graph : Information on the graphics
%
%         Output :
%                  graphics
%
% Function layout :
%
% 1. Create the arrays containing the Cartesian coordinates of the B-Spline surface grid points, the supports and the arrow vectors
%
% 2. Plot the B-Spline surfaces
%
% 3. Plot the knots on the B-Spline surfaces
%
% 4. Plot the control polygon of the B-Spline surfaces
%
% 5. Plot the supports on the geometries
%
% 6. Plot the load arrows on the geometries
%
%% Function main body

%% 0. Read input

% Initialize flags
isUndeformed = 0;
isDeformed = 0;

% Compute the control point coordinates for the deformed configuration
CPd = computeDisplacedControlPointsForIGAPlateInMembraneAction(CP,dHat);

%% 1. Create the arrays containing the Cartesian coordinates of the B-Spline surface grid points, the supports and the arrow vectors
if strcmp(graph.postprocConfig,'reference')||strcmp(graph.postprocConfig,'referenceCurrent')
    % Create the coordinates of the grid points on the surface
    [XpRef,YpRef,ZpRef] = createBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,xiGrid,etaGrid);

    % Create the coordinates of the vertices of the support triangles
    [xsRef,ysRef,zsRef] = createSupports2D(CP,rb);

    % Create the start and end points of the arrows representing the loads
    [xfRef,yfRef,zfRef] = createForceArrows2D(CP,Fl);
end
if strcmp(graph.postprocConfig,'current')||strcmp(graph.postprocConfig,'referenceCurrent')
    % Create the coordinates of the grid points on the surface
    [XpCur,YpCur,ZpCur] = createBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CPd,isNURBS,xiGrid,etaGrid);
end

%% 2. Plot the B-Spline surfaces
if strcmp(graph.postprocConfig,'reference')
    surf(XpRef,YpRef,ZpRef,'FaceColor','green','EdgeColor','none');
    hold on;
elseif strcmp(graph.postprocConfig,'current')
    surf(XpCur,YpCur,ZpCur,'FaceColor','green','EdgeColor','none');
    hold on;
elseif strcmp(graph.postprocConfig,'referenceCurrent')
    surf(XpRef,YpRef,ZpRef,'FaceColor','none','EdgeColor','none');
    hold on;
    surf(XpCur,YpCur,ZpCur,'FaceColor','green','EdgeColor','none');
end

%% 3. Plot the knots on the B-Spline surfaces
if strcmp(graph.postprocConfig,'reference')
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isUndeformed,xiGrid,etaGrid);
elseif strcmp(graph.postprocConfig,'current')
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CPd,isNURBS,isUndeformed,xiGrid,etaGrid);
elseif strcmp(graph.postprocConfig,'referenceCurrent')
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CP,isNURBS,isUndeformed,xiGrid,etaGrid);
    plot_knotsForBSplineSurfaceOnCartesianSpace(p,q,Xi,Eta,CPd,isNURBS,isDeformed,xiGrid,etaGrid);
end

%% 4. Plot the control polygon of the B-Spline surfaces
if strcmp(graph.postprocConfig,'reference')
%     plot_ControlPolygonBSplineSurface(CP);
elseif strcmp(graph.postprocConfig,'current')
    plot_ControlPolygonBSplineSurface(CPd);
elseif strcmp(graph.postprocConfig,'referenceCurrent')
%     plot_ControlPolygonBSplineSurface(CP);
    plot_ControlPolygonBSplineSurface(CPd);
end

%% 5. Plot the supports on the geometries
if strcmp(graph.postprocConfig,'reference')||strcmp(graph.postprocConfig,'referenceCurrent')
    for xiCounter =1:length(xsRef(:,1))
        plot3(xsRef(xiCounter,:),ysRef(xiCounter,:),zsRef(xiCounter,:),'Linewidth',2,'Color','black');
    end
end

%% 6. Plot the load arrows on the geometries
if strcmp(graph.postprocConfig,'reference')||strcmp(graph.postprocConfig,'referenceCurrent')
    for xiCounter =1:length(xfRef(:,1))
        plot3(xfRef(xiCounter,:),yfRef(xiCounter,:),zfRef(xiCounter,:),'Linewidth',5);
        plot3(xfRef(xiCounter,1),yfRef(xiCounter,1),zfRef(xiCounter,1),'Marker','d','MarkerFaceColor','blue','MarkerSize',10);
    end
end

end

