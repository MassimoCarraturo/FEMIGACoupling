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
function index = plot_residual_stresses_plate_with_a_hole(patch,Radius,T,d,component,graph)
%% Function documentation
%
% Plots the magnitude of the residual stress component along the area of 
% the given patch. The benchmark is the infinite plate with a circular hole 
% subject to load at minus infinity:
%
%                               exact traction
%                     -----------------------------------------------
%                     |                                             |
%                   e |                                             |
%                   x |                                             |
%                   a |                                             |
%                   c |                                             |
%                   t |                                             |
%                     |                                             |
%  <-- T_{infty}    t |                                             |
%                   r |                                             |
%                   a |                                             |
%                   c |                                    Radius   |
%                   t |                                /------------|
%                   i |                               /             
%                   o |                              /               
%                   n |                             /               
%                     |                            /                 
%                     |                            |                 
%                     -----------------------------|
%       Input : 
%       patch : Structure containing all the geometrical (NURBS parameters) 
%               and technical (parameters)
%      Radius : Radius of the circular hole
%           T : The applied load at the minus infinity
%           d : the displacement field 
%   component : The stress tensor component:
%               1 : p_11
%               2 : p_22
%               3 : p_12
%    graph : structure containing all information on the plots
%
%   Output : 
%    index : Index of the current graph
%
% Function layout :
%
% 0. Read input
%
% 1. Create the element freedom tables
%
% 2. Loop over all the sampling points
%
%    2i. Find the physical coordinates of the point on the surface 
%
%   2ii. Compute the numerical stress field
%
%  2iii. Compute the exact stress field
%
%   2iv. Compute the residual stress on the current point
%
% 3. Plot the results for the single patch
%
%% Function main body

%% 0. Read input

% Re-assign all information from the structures
p = patch.p;
q = patch.q;
U = patch.U;
V = patch.V;
parameters = patch.parameters;
CP = patch.CP;

% Read input
mu = length(U);
mv = length(V);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));
check_input(p,mu,nu,q,mv,nv);

% The material matrix
D = parameters.E/(1-parameters.nu^2)*[1 parameters.nu 0; parameters.nu 1 0; 0 0 (1-parameters.nu)/2];

%% 1. Create the element freedom tables

% dof array assigns two DoF (x,y) to every CP
% numbering follows CP: CP1->dof1,dof2 CP2->dof3,dof4
% Initialize the array of the degrees of freedome
dof = zeros(nu,nv,2);

% Initialize counter
k=1;

% Loop over all the Control points to assign global numbering to the
% degrees of freedom
for cpj = 1:nv
    for cpi = 1:nu
        dof(cpi,cpj,1)=k;
        dof(cpi,cpj,2)=k+1;
        k=k+2;
    end
end

% Assign the element displacement vectors
d_el = zeros(mu-p-1,mv-q-1,2*(p+1)*(q+1));
for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
        k=1; 
        for c = j-q-1:j-1 
            for b = i-p:i
                d_el(i,j,k)   = d(2*(c*nu+b)-1);
                d_el(i,j,k+1) = d(2*(c*nu+b));
                k=k+2;
            end
        end
    end
end

% Assign a tolerance value
tol = 10e-10;

% counting index of lines
l=1;  

%incremental step for v
r = (V(mv)-V(1))/99;

% Starting values for v
v = V(1);

%% 2. Loop over all the sampling points

while v <= V(mv)+tol
    % Find span in v-direction
    j = findspan(v,V,nv);
    
    %incremental step for u
    s = (U(mu)-U(1))/99;  
    
    % Starting values for v
    u = U(1);
    k = 1;
    while u <= U(mu)+tol
        % Find span in u-direction
        i = findspan(u,U,nu);
        
        %% 2i. Find the physical coordinates of the point on the surface 
        P(k,l,1:3) = point_on_surface(p,i,u,U,q,j,v,V,CP);
        
        % Find the actual displacement vector for the current knot span
        d_actual(:,1) = d_el(i,j,:);
        
        %% 2ii. Compute the numerical stress field
        
        % Compute the strain resultant component on the parametric coordinates
        epsilon = get_strain_vector(i,p,u,U,j,q,v,V,CP,d_actual);
        
        % Compute the complete stress vector
        sigma = D*epsilon;
        
        % Get the component
        stress_h = sigma(component);
        
        %% 2iii. Compute the exact stress field
        
        % Compute the coordinates of the point u-v on the physical space:
        % Initialize the coordinates
        x = P(k,l,1);
        y = P(k,l,2);
        
        % Compute the coordinates o the curvilinear system:
        radius = sqrt(x^2+y^2);
        theta = 2*pi()-atan(y/(-x));
        
        % preliminary values for all stress resultants:
        a=2*theta;  rr = Radius/radius;  t = T/2;
        
        % the stress resultants in the curvilinear coordinate system:
        p_r = t*(1-rr^2) + t*(1-4*rr^2+3*rr^4)*cos(a);
        p_t = t*(1+rr^2) - t*(1+3*rr^4)*cos(a);
        p_rt = -t*(1+2*rr^2-3*rr^4)*sin(a);
        
        % the stress resultants in the Cartesian coordinate system:
        if component==1
            stress = p_r*(cos(theta))^2+p_t*(sin(theta))^2-2*p_rt*sin(theta)*cos(theta);
        elseif component==2
            stress = p_r*(sin(theta))^2+p_t*(cos(theta))^2+2*p_rt*sin(theta)*cos(theta);
        elseif component==3
            stress = p_r*sin(theta)*cos(theta)-p_t*sin(theta)*cos(theta)-p_rt*((sin(theta))^2-(cos(theta))^2);
        end

        %% 2iv. Compute the residual displacement on the current point
        stress_residual(k,l,1) = norm(stress-stress_h);
        
        k=k+1;
        u=u+s;
    end
    l=l+1;
    v=v+r;
end

%% 3. Plot the results for the single patch
figure(graph.index)

% Create the colored surface
surf(P(:,:,1),P(:,:,2),P(:,:,3),stress_residual(:,:,1));

% Graphics options
shading interp;
colormap('default');

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

% On the color bar
colorbar;
hold on;

% Create edges for the patch
create_edges(p,q,U,V,CP,0,50,50);

% Adjust plotting parameters
%camlight left; 
lighting phong;
view(2);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if component==1
    title ('The magnitude of the residual stress component p_{11}');
elseif component==2
    title ('The magnitude of the residual stress component p_{22}');
elseif component==3
    title ('The magnitude of the residual stress component p_{12}');
end

hold off;

% Update the graph index
index = graph.index + 1;

end

