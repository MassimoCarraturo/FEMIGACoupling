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
function  plot_exactResultantRectangularPlateWithAHole(patch,T,Radius,npoints,graph,clr)
%% Function documentation
%
% Plots the exact resultant (displacement,strain,stress) for a rectangular
% plate with a hole for which analytical solution is known
%
%        Input : 
%        patch : structure containing all the technical, geometrical and
%                displacement information of the membrane
%            l : the size of the pressure loaded edge
%            h : the size of the fixed edges
%            T : the magnitude of the load at x = - inf
%       Radius : The radius of the quarter of a circle shaped hole
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
G2orth = zeros(3,1);

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
    spanv = findspan(uv,V,nv);
    
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
        x = curve(i,1);
        y = curve(i,2);
        
        % Compute the coordinates o the curvilinear system:
        r = sqrt(x^2+y^2);
        theta = 2*pi()-atan(y/(-x));
        
        % Shift the point to the origin
        % curve_origin(i,1:3) = curve(i,1:3) - [0 curve(i,2) 0];
        curve_origin(i,1:3) = curve(i,1:3);
        
        % Compute the exact resultants on the point of the curve
        if graph.resultant==1
            % Parameters:
            mue = E/(1+nue)/2;
            kappa = (3-nue)/(1+nue);
            
            % Choose displacement component
            if graph.component==1
                resultant(i,1) = -T*Radius*(r*(kappa+1)*cos(theta)/Radius+2*Radius*((1+kappa)*cos(theta)+cos(3*theta))/r-2*Radius^3*cos(3*theta)/r^3)/mue/8;
            elseif graph.component==2
                resultant(i,1) = -T*Radius*(r*(kappa-3)*sin(theta)/Radius+2*Radius*((1-kappa)*sin(theta)+sin(3*theta))/r-2*Radius^3*sin(3*theta)/r^3)/mue/8;
            elseif graph.component==3
                displ = zeros(2,1);
                displ(1,1) = -T*Radius*(r*(kappa+1)*cos(theta)/Radius+2*Radius*((1+kappa)*cos(theta)+cos(3*theta))/r-2*Radius^3*cos(3*theta)/r^3)/mue/8;
                displ(2,1) = -T*Radius*(r*(kappa-3)*sin(theta)/Radius+2*Radius*((1-kappa)*sin(theta)+sin(3*theta))/r-2*Radius^3*sin(3*theta)/r^3)/mue/8;
                resultant(i,1) = norm(displ);
            end
        elseif graph.resultant==2
            % preliminary values for all stress resultants:
            a=2*theta;  rr = Radius/r;  t = T/2;
            % the stress resultants in the curvilinear coordinate system:
            p_r = t*(1-rr^2) + t*(1-4*rr^2+3*rr^4)*cos(a);
            p_t = t*(1+rr^2) - t*(1+3*rr^4)*cos(a);
            p_rt = -t*(1+2*rr^2-3*rr^4)*sin(a);
            % the stress resultants in the Cartesian coordinate system:
            p_x = p_r*(cos(theta))^2+p_t*(sin(theta))^2-2*p_rt*sin(theta)*cos(theta);
            p_y = p_r*(sin(theta))^2+p_t*(cos(theta))^2+2*p_rt*sin(theta)*cos(theta);
            p_xy = p_r*sin(theta)*cos(theta)-p_t*sin(theta)*cos(theta)-p_rt*((sin(theta))^2-(cos(theta))^2);
            
            % Choose strain component
            if graph.component==1
                resultant(i,1) = (p_x-nue*p_y)/E/th;
           elseif graph.component==2
                resultant(i,1) = (p_y-nue*p_x)/E/th;
            elseif graph.component==3
                % The shear modulus:
                G = E/(1-nue)/2;
                resultant(i,1) = p_xy/2/G/th;
            end
        elseif graph.resultant==3
            % preliminary values for all stress resultants:
            a=2*theta;  rr = Radius/r;  t = T/2;
            % the stress resultants in the curvilinear coordinate system:
            p_r = t*(1-rr^2) + t*(1-4*rr^2+3*rr^4)*cos(a);
            p_t = t*(1+rr^2) - t*(1+3*rr^4)*cos(a);
            p_rt = -t*(1+2*rr^2-3*rr^4)*sin(a);
            
            % choose the stress component:
            if graph.component==1
                resultant(i,1) = p_r*(cos(theta))^2+p_t*(sin(theta))^2-2*p_rt*sin(theta)*cos(theta);
            elseif graph.component==2
                resultant(i,1) = p_r*(sin(theta))^2+p_t*(cos(theta))^2+2*p_rt*sin(theta)*cos(theta);
            elseif graph.component==3
                resultant(i,1) = p_r*sin(theta)*cos(theta)-p_t*sin(theta)*cos(theta)-p_rt*((sin(theta))^2-(cos(theta))^2);
            end
        end
        
        % Compute base vector in the u-parametric coordinates 
        [G1,~] = base_vectors2D(p,spanu,u,U,q,spanv,uv,V,CP);
        
        % Compute an orthogonal to G2 vector
        G1orth(1,1) = -G1(2,1);
        G1orth(2,1) = G1(1,1);
        G1orth(3,1) = G1(3,1);
        
        % Normalize the vector to give physical perspective
        t = G1orth/norm(G1orth);
        
        test1 = curve_origin(i,:)';
        test2 = resultant(i,1);
        % Orthogonally projected from the curve the physical strain point
        P(i,:) = test1 + t*test2;
        
        % Update coordinate
        u = u + du;
    else
        % Find the knot span where we are inside
        spanv = findspan(v,V,nv);
        
        % Compute the coordinates of the point on curve
        curve(i,1:3) = point_on_surface(p,spanu,uv,U,q,spanv,v,V,CP);
        x = curve(i,1);
        y = curve(i,2);
        
        % Compute the coordinates o the curvilinear system:
        r = sqrt(x^2+y^2);
        theta = 2*pi()-atan(y/(-x));
        
        % Shift the point to the origin
        % curve_origin(i,1:3) = curve(i,1:3) - [0 curve(i,2) 0];
        curve_origin(i,1:3) = curve(i,1:3);
        
        % Map the coordinates of the point in the old coordinate system
        % into the new coordinate system
        % curve_new = curve(i,1:3)-Rx;
        
        % Compute the exact resultants on the point of the curve
        if graph.resultant==1
            % Parameters:
            mue = E/(1+nue)/2;
            kappa = (3-nue)/(1+nue);
            
            % Choose displacement component
            if graph.component==1
                resultant(i,1) = -T*Radius*(r*(kappa+1)*cos(theta)/Radius+2*Radius*((1+kappa)*cos(theta)+cos(3*theta))/r-2*Radius^3*cos(3*theta)/r^3)/mue/8;
            elseif graph.component==2
                resultant(i,1) = -T*Radius*(r*(kappa-3)*sin(theta)/Radius+2*Radius*((1-kappa)*sin(theta)+sin(3*theta))/r-2*Radius^3*sin(3*theta)/r^3)/mue/8;
            elseif graph.component==3
                displ = zeros(2,1);
                displ(1,1) = -T*Radius*(r*(kappa+1)*cos(theta)/Radius+2*Radius*((1+kappa)*cos(theta)+cos(3*theta))/r-2*Radius^3*cos(3*theta)/r^3)/mue/8;
                displ(2,1) = -T*Radius*(r*(kappa-3)*sin(theta)/Radius+2*Radius*((1-kappa)*sin(theta)+sin(3*theta))/r-2*Radius^3*sin(3*theta)/r^3)/mue/8;
                resultant(i,1) = norm(displ);
            end
        elseif graph.resultant==2
            % preliminary values for all stress resultants:
            a=2*theta;  rr = Radius/r;  t = T/2;
            % the stress resultants in the curvilinear coordinate system:
            p_r = t*(1-rr^2) + t*(1-4*rr^2+3*rr^4)*cos(a);
            p_t = t*(1+rr^2) - t*(1+3*rr^4)*cos(a);
            p_rt = -t*(1+2*rr^2-3*rr^4)*sin(a);
            % the stress resultants in the Cartesian coordinate system:
            p_x = p_r*(cos(theta))^2+p_t*(sin(theta))^2-2*p_rt*sin(theta)*cos(theta);
            p_y = p_r*(sin(theta))^2+p_t*(cos(theta))^2+2*p_rt*sin(theta)*cos(theta);
            p_xy = p_r*sin(theta)*cos(theta)-p_t*sin(theta)*cos(theta)-p_rt*((sin(theta))^2-(cos(theta))^2);
            
            % Choose strain component
            if graph.component==1
                resultant(i,1) = (p_x-nue*p_y)/E/th;
           elseif graph.component==2
                resultant(i,1) = (p_y-nue*p_x)/E/th;
            elseif graph.component==3
                % The shear modulus:
                G = E/(1-nue)/2;
                resultant(i,1) = p_xy/2/G/th;
            end
        elseif graph.resultant==3
            % preliminary values for all stress resultants:
            a=2*theta;  rr = Radius/r;  t = T/2;
            % the stress resultants in the curvilinear coordinate system:
            p_r = t*(1-rr^2) + t*(1-4*rr^2+3*rr^4)*cos(a);
            p_t = t*(1+rr^2) - t*(1+3*rr^4)*cos(a);
            p_rt = -t*(1+2*rr^2-3*rr^4)*sin(a);
            
            % choose the stress component:
            if graph.component==1
                resultant(i,1) = p_r*(cos(theta))^2+p_t*(sin(theta))^2-2*p_rt*sin(theta)*cos(theta);
            elseif graph.component==2
                resultant(i,1) = p_r*(sin(theta))^2+p_t*(cos(theta))^2+2*p_rt*sin(theta)*cos(theta);
            elseif graph.component==3
                resultant(i,1) = p_r*sin(theta)*cos(theta)-p_t*sin(theta)*cos(theta)-p_rt*((sin(theta))^2-(cos(theta))^2);
            end
        end
        
        % Compute base vector in the u-parametric coordinates 
        [~,G2] = base_vectors2D(p,spanu,uv,U,q,spanv,v,V,CP);
        
        % Compute an orthogonal to G2 vector
        G2orth(1,1) = G2(2,1);
        G2orth(2,1) = -G2(1,1);
        G2orth(3,1) = G2(3,1);
        
        % Normalize the vector to give physical perspective
        t = G2orth/norm(G2orth);
        
        % Orthogonally projected from the curve the physical strain point
        P(i,:) = curve_origin(i,:)' + t*resultant(i,1);
        
        % Update coordinate
        v = v + dv;
    end
end

figure(graph.index)

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