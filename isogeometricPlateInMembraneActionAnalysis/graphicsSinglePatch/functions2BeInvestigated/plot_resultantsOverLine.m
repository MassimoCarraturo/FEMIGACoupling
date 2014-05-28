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
function plot_resultantsOverLine(patch,npoints,graph,clr)
%% Function documentation
%
% Given the parametric line region and the desirable component of the
% strain plots the strain distribution over the line region in the physical
% space
%
%        Input : 
%        patch : structure containing all the technical, geometrical and
%                displacement information of the membrane
%      npoints : number of sampling points
%        graph : structure containing all information on the plots
%
%       Output : Graphics
%
%    Function layout :
%
% 1. Loop over all the sampling points of the coupling curve
% 2. Compute the base vector corresponding to the u- or v-parametric curve
%    at each sampling point 
% 3. Compute the strain vector at each sampling points and extract the
%    chosen component
% 4. Project the value orthogonal to the parametric line
% 5. Plot the parametric line together with the strain values
%
%
%% Function main body

% Read NURBS and coupling information
p = patch.p;
U = patch.U;
q = patch.q;
V = patch.V;
CP = patch.CP;
parameters = patch.parameters;
ucoup = patch.ucoup;
vcoup = patch.vcoup;

% Initialize auxiliary arrays
curve = zeros(npoints+1,1);
curve_origin = zeros(npoints+1,1);
resultant = zeros(npoints+1,1);
P = zeros(npoints+1,3);
Gorth = zeros(3,1);

% Compute the plane stress material matrix needed for the resultants
D = parameters.E/(1-parameters.nu^2)*[1 parameters.nu 0; parameters.nu 1 0; 0 0 (1-parameters.nu)/2];

% Get the displacement vector that corresponds to the parametric line
d = patch.d;

% Read input
mu = length(U);
mv = length(V);
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));
check_input(p,mu,nu,q,mv,nv);

% Distribute the displacement vectors over the elements
d_el = zeros(mu-p-1,mv-q-1,2*(p+1)*(q+1));
for j = (q+1):(mv-q-1)
    for i = (p+1):(mu-p-1)
        % counter
        k=1; 
        for c = j-q-1:j-1 
            for b = i-p:i
                d_el(i,j,k)   = d(2*(c*nu+b)-1);
                d_el(i,j,k+1) = d(2*(c*nu+b));
                % update counter
                k=k+2;
            end
        end
    end
end

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
        % curve_origin(i,1:3) = curve(i,1:3) - [curve(i,1) 0 0];
        curve_origin(i,1:3) = curve(i,1:3);
        
        % The displacement vector for the current element knot span
        d_actual(:,1) = d_el(spanu,spanv,:);
                
        % Compute the displacement field on the parametric coordinates
        displacement = get_displacement_vector(spanu,p,u,U,spanv,q,uv,V,CP,d_actual);
        
        % Compute the resultant component on the parametric coordinates
        epsilon = get_strain_vector(spanu,p,u,U,spanv,q,uv,V,CP,d_actual);
        sigma = D*epsilon;
        
        if graph.resultant==1
            if graph.component==1 || graph.component==2
                resultant(i,1) = displacement(graph.component);
            elseif graph.component==3
                resultant(i,1) = (displacement(1)^2 + displacement(2)^2)^(1/2);
            end
        elseif graph.resultant==2
            if graph.component~=4 && graph.component~=5
                resultant(i,1) = epsilon(graph.component);
            elseif graph.component==4
                resultant(i,1) = 0.5*(epsilon(1,1)+epsilon(2,1) + sqrt((epsilon(1,1)-epsilon(2,1))^2 + epsilon(3,1)^2));
            elseif graph.component==5
                resultant(i,1) = 0.5*(epsilon(1,1)+epsilon(2,1) - sqrt((epsilon(1,1)-epsilon(2,1))^2 + epsilon(3,1)^2));
            end
        elseif graph.resultant==3
            if graph.component~=4 && graph.component~=5
                resultant(i,1) = sigma(graph.component);
            elseif graph.component==4
                resultant(i,1) = 0.5*(sigma(1,1)+sigma(2,1) + sqrt((sigma(1,1)-sigma(2,1))^2 + 4*sigma(3,1)^2));
            elseif graph.component==5
                resultant(i,1) = 0.5*(sigma(1,1)+sigma(2,1) - sqrt((sigma(2,1)-sigma(2,1))^2 + 4*sigma(3,1)^2));
            end
        end
        
        % Compute base vector in the u-parametric coordinates 
        [G1,~] = base_vectors2D(p,spanu,u,U,q,spanv,uv,V,CP);
        
        % Compute an orthogonal to G1 vector
        Gorth(1) = G1(2);
        Gorth(2) = -G1(1);
        Gorth(3) = G1(3);
        
        % Normalize the vector to give physical perspective
        t = -Gorth/norm(Gorth);
        
        % Orthogonally projected from the curve the physical strain point
        P(i,:) = curve_origin(i,:)' + t*resultant(i,1);
        
        % Update coordinate
        u = u + du;
    else
        % Find the knot span where we are inside
        spanv = findspan(v,V,nv);
        
        % Compute the point on the curve
        curve(i,1:3) = point_on_surface(p,spanu,uv,U,q,spanv,v,V,CP);
        
        % Shift the point to the origin
        % curve_origin(i,1:3) = curve(i,1:3) - [curve(i,1) 0 0];
        curve_origin(i,1:3) = curve(i,1:3);
        
        % The displacement vector for the current element knot span
        d_actual(:,1) = d_el(spanu,spanv,:);
           
        % Compute the displacement field on the parametric coordinates
        displacement = get_displacement_vector(spanu,p,uv,U,spanv,q,v,V,CP,d_actual);
        
        % Compute the resultant component on the parametric coordinates
        epsilon = get_strain_vector(spanu,p,uv,U,spanv,q,v,V,CP,d_actual);
        sigma = D*epsilon;
        
        if graph.resultant==1
            if graph.component==1 || graph.component==2
                resultant(i,1) = displacement(graph.component);
            elseif graph.component==3
                resultant(i,1) = (displacement(1)^2 + displacement(2)^2)^(1/2);
            end
        elseif graph.resultant==2
            if graph.component~=4 && graph.component~=5
                resultant(i,1) = epsilon(graph.component);
            elseif graph.component==4
                resultant(i,1) = 0.5*(epsilon(1,1)+epsilon(2,1) + sqrt((epsilon(1,1)-epsilon(2,1))^2 + epsilon(3,1)^2));
            elseif graph.component==5
                resultant(i,1) = 0.5*(epsilon(1,1)+epsilon(2,1) - sqrt((epsilon(1,1)-epsilon(2,1))^2 + epsilon(3,1)^2));
            end
        elseif graph.resultant==3
            if graph.component~=4 && graph.component~=5
                resultant(i,1) = sigma(graph.component);
            elseif graph.component==4
                resultant(i,1) = 0.5*(sigma(1,1)+sigma(2,1) + sqrt((sigma(1,1)-sigma(2,1))^2 + 4*sigma(3,1)^2));
            elseif graph.component==5
                resultant(i,1) = 0.5*(sigma(1,1)+sigma(2,1) - sqrt((sigma(2,1)-sigma(2,1))^2 + 4*sigma(3,1)^2));
            end
        end
        
        % Compute base vector in the u-parametric coordinates 
        [~,G2] = base_vectors2D(p,spanu,uv,U,q,spanv,v,V,CP);
        
        % Compute an orthogonal to G2 vector
        Gorth(1,1) = G2(2);
        Gorth(2,1) = -G2(1);
        Gorth(3,1) = G2(3);
        
        % Normalize the vector to give physical perspective
        t = Gorth/norm(Gorth);
        
        % Orthogonally projected from the curve the physical strain point
        P(i,:) = curve_origin(i,:)' + t*resultant(i,1);
        
        % Update coordinate
        v = v + dv;
    end
end

if graph.resultant==1 || graph.resultant==2 
    % Plot the coupling interface curve at the origin
    plot3(curve_origin(:,1),curve_origin(:,2),curve_origin(:,3),'LineWidth',2,'Color','green');

    hold on;

    % Plot the resultant distribution over the coupling interface 
    plot3(P(:,1),P(:,2),P(:,3),'Color',clr);
elseif graph.resultant==3
    % Plot the resultant distribution over the coupling interface 
    plot3(P(:,1),P(:,2),P(:,3),'Color',clr);

    hold on;
    
    % Plot the coupling interface curve at the origin
    plot3(curve_origin(:,1),curve_origin(:,2),curve_origin(:,3),'LineWidth',2,'Color','green');
end

camlight left; lighting phong;
view(2);
grid on;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
title('Distribution over the coupling interface');

end