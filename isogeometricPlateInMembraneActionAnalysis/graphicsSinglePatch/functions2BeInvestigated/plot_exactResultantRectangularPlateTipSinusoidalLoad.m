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
function plot_exact_resultant_rectangular_membrane_tip_sinusoidal(membrane,l,h,ql,npoints,graph,clr)
%% Function documentation
%
% Plots the exact resultant (displacement,strain,stress) for a rectangular
% plate fixed at one edge, subject to tip sinusoidal load at the free edge
%
%        Input : 
%     membrane : structure containing all the technical, geometrical and
%                displacement information of the membrane
%            l : the size of the pressure loaded edge
%            h : the size of the fixed edges
%           ql : the magnitude of the pressure load (+ downwards - upwards)
%      npoints : number of sampling points
%        graph : structure containing all information on the plots
%
%       Output : Graphics
%
%    Function layout :
%
%
%
%% Function main body

% Read NURBS and coupling information
p = membrane.p;
U = membrane.U;
q = membrane.q;
V = membrane.V;
CP = membrane.CP;
parameters = membrane.parameters;
ucoup = membrane.ucoup;
vcoup = membrane.vcoup;

% Initialize auxiliary arrays
curve = zeros(npoints+1,1);
curve_origin = zeros(npoints+1,1);
resultant = zeros(npoints+1,1);
P = zeros(npoints+1,3);
G2orth = zeros(2,1);

% Assign the plane stress material matrix needed for the resultants
E = parameters.E;
th = parameters.t;
nue = parameters.nu;

%D = parameters.E/(1-parameters.nu^2)*[1 parameters.nu 0; parameters.nu 1 0; 0 0 (1-parameters.nu)/2];

% Read input
mu = length(U);
mv = length(V);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));
check_input(p,mu,nu,q,mv,nv);

if vcoup(1)==vcoup(2)
    % The coupling line is on a u-curve
    uv = vcoup(1);
    
    % The span in v-direction is fixed
    spanv = findspan(uv,U,nv);
    
    % Initialize parametric coordinate
    u = ucoup(1);
    
    % Compute the incremental step
    du = (ucoup(2)-ucoup(1))/npoints;
else
    % The coupling line is on a v-curve
    uv = ucoup(1);
    
    % The span in u-direction is fixed
    spanu = findspan(uv,U,nu);
    
    % Initialize parametric coordinate
    v = vcoup(1);
    
    % Compute the incremental step
    dv = (vcoup(2)-vcoup(1))/npoints;
end

% Loop over all the sampling points of the coupling curve
for i=1:npoints+1
    if vcoup(1)==vcoup(2)
        % Find the knot span where we are inside
        spanu = findspan(u,U,nu);
        
        % Compute the point on the curve
        curve(i,1:3) = point_on_surface(p,spanu,u,U,q,spanv,uv,V,CP);
        
        % Shift the point to the origin
        curve_origin(i,1:3) = curve(i,1:3) - [curve(i,1) 0 0];
        
        % shift the coordinates to the center of the plate
        % The position vector of the new coordinate system with respect to
        % the old coordinate system:
        Rx = [0 0 0];
        
        % Map the coordinates of the point in the old coordinate system
        % into the new coordinate system
        curve_new = curve(i,1:3)-Rx;
        
        % Compute the exact resultants on the point of the curve
        if graph.resultant==1
            if graph.component==1
                resultant(i,1) = 0;
            elseif graph.component==2
                resultant(i,1) = 0;
            end
        elseif graph.resultant==2
            if graph.component==1
                resultant(i,1) = 0;
           elseif graph.component==2
                resultant(i,1) = 0;
            elseif graph.component==3
                resultant(i,1) = 0;
            end
        elseif graph.resultant==3
            if graph.component==1
                resultant(i,1) = ql*sin(pi*curve_new(2)/h)*cos(pi*curve_new(2)/h);
            elseif graph.component==2
                resultant(i,1) = 0;
            elseif graph.component==3
                resultant(i,1) = 0;
            end
        end
        
        % Compute base vector in the u-parametric coordinates 
        [G1,~] = base_vectors2D(p,spanu,u,U,q,spanv,uv,V,CP);
        
        % Compute an orthogonal to G2 vector
        G1orth(1) = G1(2);
        G1orth(2) = -G1(1);
        G1orth(3) = G1(3);
        
        % Normalize the vector to give physical perspective
        t = G1orth/norm(G1orth);
        
        % Orthogonally projected from the curve the physical strain point
        P(i,:) = curve_origin(i,:)' + t*resultant(i,1);
        
        % Update coordinate
        u = u + du;
    else
        % Find the knot span where we are inside
        spanv = findspan(v,V,nv);
        
        % Compute the coordinates of the point on curve
        curve(i,1:3) = point_on_surface(p,spanu,uv,U,q,spanv,v,V,CP);
           
        % Shift the point to the origin
        curve_origin(i,1:3) = curve(i,1:3) - [curve(i,1) 0 0];
        
        % shift the coordinates to the center of the plate
        % The position vector of the new coordinate system with respect to
        % the old coordinate system:
        Rx = [0 0 0];
        
        % Map the coordinates of the point in the old coordinate system
        % into the new coordinate system
        curve_new = curve(i,1:3)-Rx;
        
        % Compute the exact resultants on the point of the curve
        if graph.resultant==1
            if graph.component==1
                resultant(i,1) = 0;
            elseif graph.component==2
                resultant(i,1) = 0;
            end
        elseif graph.resultant==2
            if graph.component==1
                resultant(i,1) = 0;
           elseif graph.component==2
                resultant(i,1) = 0;
            elseif graph.component==3
                resultant(i,1) = 0;
            end
        elseif graph.resultant==3
            if graph.component==1
                resultant(i,1) = ql*sin(pi*curve_new(2)/h)*cos(pi*curve_new(2)/h);
            elseif graph.component==2
                resultant(i,1) = 0;
            elseif graph.component==3
                resultant(i,1) = 0;
            end
        end
        
        % Compute base vector in the u-parametric coordinates 
        [~,G2] = base_vectors2D(p,spanu,uv,U,q,spanv,v,V,CP);
        
        % Compute an orthogonal to G2 vector
        G2orth(1) = G2(2);
        G2orth(2) = -G2(1);
        G2orth(3) = G2(3);
        
        % Normalize the vector to give physical perspective
        t = G2orth/norm(G2orth);
        
        % Orthogonally projected from the curve the physical strain point
        P(i,:) = curve_origin(i,:)' + t*resultant(i,1);
        
        % Update coordinate
        v = v + dv;
    end
end

if graph.resultant==1 || graph.resultant==2 
    % Plot the coupling interface curve
    plot3(curve_origin(:,1),curve_origin(:,2),curve_origin(:,3),'LineWidth',2,'Color','green');

    hold on;

    % Plot the resultant distribution over the coupling interface 
    plot3(P(:,1),P(:,2),P(:,3),'Color',clr);
elseif graph.resultant==3
    % Plot the resultant distribution over the coupling interface 
    plot3(P(:,1),P(:,2),P(:,3),'Color',clr);

    hold on;
    
    % Plot the coupling interface curve
    plot3(curve_origin(:,1),curve_origin(:,2),curve_origin(:,3),'LineWidth',2,'Color','green');
end

camlight left; lighting phong;
view(2);
grid on;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title('Distribution over the coupling interface');

end