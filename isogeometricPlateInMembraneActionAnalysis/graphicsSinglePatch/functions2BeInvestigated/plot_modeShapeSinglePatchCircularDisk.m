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
function index = plot_modeShapeSinglePatchCircularDisk(p,q,U,V,CP,parameters,eigenFrequencies,modeShapes,modeNumber,graph)
%% Function documentation
%
% Plots the mode shapes for one plane stress NURBS patch together
% with the resulting stress field over the plate's domain. The stress
% distribution field is analogous to the scaling of the eigenVectors which
% means that the distribution is not physical but dependent on the choice
% of scaling of the eigeVectors. This function is a specialization of the
% general function plot_mode_shape_single_patch.m handling in addition the
% parametrization of a circular disk modelled with NURBS for which there
% exist elements with zero Jacobian determinant
%
%                Input : 
%                  p,q : The polynomial degrees of the NURBS surface
%                  U,V : The knot vectors in U,V-direction
%                   CP : The set of control points and weights
%     eigenFrequencies : The vector of the eigenfrequencies
%           parameters : Technical and geometrical parameters of the plate
%           modeShapes : The eigenvectors/mode shapes of the patch
%           modeNumber : The number of the mode shape to be plotted
%                graph : structure containing all information on the plots
%
%               Output :
%                index : the index of the graph
%
%% Function main body

% Compute the control point coordinates for the deformed configuration
CPd = displaced_control_points(CP,modeShapes(:,modeNumber));

% Number of knots in u,v-direction
mu = length(U);
mv = length(V);

% Number of Control Points in u,v-direction
nu = length(CPd(:,1,1));
nv = length(CPd(1,:,1));

% Material matrix needed to compute the resultants
D = parameters.E/(1-parameters.nu^2)*[1 parameters.nu 0; parameters.nu 1 0; 0 0 (1-parameters.nu)/2];

% Compute point coordinates for the original geometry and the supports
if graph.undeformed==1
    [Xp,Yp,Zp] = create_surface(p,q,U,V,CP,50,50);
end
if graph.deformed==1
    [Xp_def,Yp_def,Zp_def] = create_surface(p,q,U,V,CPd,50,50);
end

% Number the current plot
figure(graph.index)

% FIRST WINDOW: UNDEFORMED/DEFORMED
subplot(2,1,1);

% graph counter
access = 0;

% Plot the initial geometry
if graph.undeformed==1 && graph.deformed==0

    surf(Xp,Yp,Zp,'FaceColor','green','EdgeColor','none');
    if access == 0
        hold;
    end
     
    access = access + 1;

    % Plot the element edges for the patch
    create_edges(p,q,U,V,CP,0,50,50)
    
    % control points and control polygon for the patch
    create_control_polygon(CP)
    
    % Plot the supports for the patch
    for k = 1:length(xs(:,1))
        plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
    end
elseif graph.undeformed==1 && graph.deformed==1
    surf(Xp,Yp,Zp,'FaceColor','green','EdgeColor','none','FaceAlpha',0.5);
    if access == 0
        hold;
    end
     
    access = access + 1;
    % Plot the element edges for patch
    create_edges(p,q,U,V,CP,1,50,50)
end

% Plot the deformed geometry
if graph.deformed==1
    if access == 0
        surf(Xp_def,Yp_def,Zp_def,'FaceColor','green','EdgeColor','none');
        hold;
        
        % Draw the knots
        create_edges(p,q,U,V,CPd,0,50,50)
        
        % Draw the control polygon
        create_control_polygon(CPd)
        
        axis equal;
    else
        % Plot the NURBS surface
        surf(Xp_def,Yp_def,Zp_def,'FaceColor','green','EdgeColor','none');
        
        % Plot the knots on the surface
        create_edges(p,q,U,V,CPd,1,50,50)
    end
end

% Adjust plotting parameters
camlight left; 
lighting phong;
view(2);
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
% title ('Reference and/or current geometry modeShape No. = %d and eigenFrequency = %d',modeNumber,eigenFrequencies);
title ('Reference and/or current geometry');
hold off;

% SECOND WINDOW: STRAIN/STRESS ON UNDEFORMED STRUCTURE

% element displacement vectors for patch 1
d_el = zeros(mu-p-1,mv-q-1,2*(p+1)*(q+1));
for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
        k=1; 
        for c = j-q-1:j-1 
            for b = i-p:i
                d_el(i,j,k)   = modeShapes(2*(c*nu+b)-1,modeNumber);
                d_el(i,j,k+1) = modeShapes(2*(c*nu+b),modeNumber);
                
                % Update counter
                k=k+2;
            end
        end
    end
end

% Assign a tolerance value
tol = 10e-10;

% On the stress/strain resultants

% For the 1st patch

% counting index of lines
l=1;  
%incremental step for v
r = (V(mv)-V(1))/99;    
v = V(1);
while v <= V(mv)+tol
    j = findspan(v,V,nv);
    s = (U(mu)-U(1))/99;  %incremental step for u
    u = U(1);
    k = 1;
    while u <= U(mu)+tol
        i = findspan(u,U,nu);
        P(k,l,1:3) = point_on_surface(p,i,u,U,q,j,v,V,CP);
        d_actual(:,1) = d_el(i,j,:);
        
        % Compute the displacement components 
        disp_actual(k,l,1:2) = get_displacement_vector(i,p,u,U,j,q,v,V,CP,d_actual);
        
        % Compute the total magnitude of the displacement
        if graph.resultant==1 && graph.component==3
            total_disp_actual(k,l) = (disp_actual(k,l,1)^2 + disp_actual(k,l,2)^2)^(1/2);
        end
        
        % Compute NURBS basis functions
        [~,dR] = nurbs_basis_functions_and_first_derivatives_2D(i,p,u,U,j,q,v,V,CP);

        % Jacobian of transformation from the physical space to the NURBS parent
        % domain
        J = zeros(2,2);

        % Initialize counter
        kJ = 0;

        for c = 0:q
            for b = 0:p
                % Update counter
                kJ = kJ + 1;
                
                J(1,1) = J(1,1) + CP(i-p+b,j-q+c,1)*dR(kJ,1);
                J(1,2) = J(1,2) + CP(i-p+b,j-q+c,2)*dR(kJ,1);
                J(2,1) = J(2,1) + CP(i-p+b,j-q+c,1)*dR(kJ,2);
                J(2,2) = J(2,2) + CP(i-p+b,j-q+c,2)*dR(kJ,2);
            end
        end
       
        % Compute the stresses only if the element is not distorted
        if det(J)~=0
            % Compute the strain components
            eps_actual = get_strain_vector(i,p,u,U,j,q,v,V,CP,d_actual);

            eps(k,l,1:3) = eps_actual(1:3);
            eps(k,l,4) = 0.5*(eps(k,l,1)+eps(k,l,2) + sqrt((eps(k,l,1)-eps(k,l,2))^2 + eps(k,l,3)^2));
            eps(k,l,5) = 0.5*(eps(k,l,1)+eps(k,l,2) - sqrt((eps(k,l,1)-eps(k,l,2))^2 + eps(k,l,3)^2));
            sig(k,l,1:3) = D*eps_actual(1:3);
            sig(k,l,4) = 0.5*(sig(k,l,1)+sig(k,l,2) + sqrt((sig(k,l,1)-sig(k,l,2))^2 + 4*sig(k,l,3)^2));
            sig(k,l,5) = 0.5*(sig(k,l,1)+sig(k,l,2) - sqrt((sig(k,l,1)-sig(k,l,2))^2 + 4*sig(k,l,3)^2));

            % Update internal counter
            k=k+1;
        end
        
        u=u+s;
    end
    l=l+1;
    v=v+r;
end



% plot
subplot(2,1,2);

if graph.resultant==1
    if graph.component == 3
        surf(P(:,:,1),P(:,:,2),P(:,:,3),total_disp_actual(:,:),'EdgeColor','none');
    else
        surf(P(:,:,1),P(:,:,2),P(:,:,3),disp_actual(:,:,graph.component),'EdgeColor','none');
    end   
elseif graph.resultant==2
    surf(P(:,:,1),P(:,:,2),P(:,:,3),eps(:,:,graph.component),'EdgeColor','none');
elseif graph.resultant==3
    surf(P(:,:,1),P(:,:,2),P(:,:,3),sig(:,:,graph.component),'EdgeColor','none'); 
end

% Graphics options
% shading interp;
colormap('default');
view(2);
axis equal;
grid off;

% invert default colormap => red=negativ, blue=positive
% COL=colormap;
% invCOL(:,1)=COL(:,3);
% invCOL(:,2)=COL(:,2);
% invCOL(:,3)=COL(:,1);
% colormap(invCOL);
% 
% % make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);

colorbar;
hold on;

% Draw the knots on the geometry
create_edges(p,q,U,V,CP,0,50,50);

view(2);
% axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if graph.resultant==1 && graph.component==1 
    title('displacement component u_x');
elseif graph.resultant==1 && graph.component==2 
    title('displacement component u_x');
elseif graph.resultant==1 && graph.component==3 
    title('displacement magnitude u=(u_x^2+u_y^2)^{0.5}');
elseif graph.resultant==2 && graph.component==1 
    title('strain epsilon xx');
elseif graph.resultant==2 && graph.component==2 
    title('strain epsilon yy');
elseif graph.resultant==2 && graph.component==3 
    title('strain epsilon xy');
elseif graph.resultant==2 && graph.component==4 
    title('principal strain e1');
elseif graph.resultant==2 && graph.component==5 
    title('principal strain e2');
elseif graph.resultant==3 && graph.component==1
    title('stress sigma xx');
elseif graph.resultant==3 && graph.component==2
    title('stress sigma yy');
elseif graph.resultant==3 && graph.component==3 
    title('stress sigma xy');
elseif graph.resultant==3 && graph.component==4
    title('principal stress s1');
elseif graph.resultant==3 && graph.component==5 
    title('principal stress s2');
end

hold off;

% Update plot index by 1
index = graph.index+1;

end

