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
function index = plot_residualForces(patch,d,graph)
%% Function documentation
%
% Plots the x- or y-component of the the residual force e = div(p)-b
% obtained by the momentum equation (colored plot)
% 
%    Input : 
%    patch : Structure containing all the geometrical (NURBS parameters) 
%            and technical (parameters)
%        d : the displacement field obtained by the postprocessing
%    graph : structure containing all information on the plots
%
%   Output : graphics
%
%% Function main body

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

% dof array assigns two DoF (x,y) to every CP
% numbering follows CP: CP1->dof1,dof2 CP2->dof3,dof4
% Initialize the array of the degrees of freedome
dof = zeros(nu,nv,2);

% Number of control points (~basis functions) that affect each element
nNode_e = (p+1)*(q+1);

% Number of degrees of freedom per element
ndof_e = 2*nNode_e;
            
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

% Starting the loop over all the Gauss points

% Assign a tolerance value
tol = 10e-10;

% counting index of lines
l=1;  

%incremental step for v
r = (V(mv)-V(1))/99;

% Starting values for v
v = V(1);

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
        
        % Find the physical coordinates of the point on the surface 
        P(k,l,1:3) = point_on_surface(p,i,u,U,q,j,v,V,CP);
        
        % Find the actual displacement vector for the current knot span
        d_actual(:,1) = d_el(i,j,:);
        
        % Compute the NURBS basis functions, its first and 
        % second derivatives at the quadrature point
        [~,dR,ddR] = nurbs_basis_functions_first_and_second_derivatives_2D(i,p,u,U,j,q,v,V,CP);

        % Compute the Hessian matrix which is needed for the 
        % tranformation of the second derivatives

        % Initialize Jacobian matrix
        Jxxi = zeros(2,2);

        % Initialize Hessian matrix
        Hxxi = zeros(2,3);

        % initialize counter
        ki = 0;

        % Loop over all the non-zero contributions at the span
        % under study
        for c = 0:q
            for b = 0:p
                % Update counter
                ki = ki + 1;

                % Compute recursively the entries of the Jacobian
                Jxxi(1,1) = Jxxi(1,1) + CP(i-p+b,j-q+c,1)*dR(ki,1);
                Jxxi(1,2) = Jxxi(1,2) + CP(i-p+b,j-q+c,2)*dR(ki,1);
                Jxxi(2,1) = Jxxi(2,1) + CP(i-p+b,j-q+c,1)*dR(ki,2);
                Jxxi(2,2) = Jxxi(2,2) + CP(i-p+b,j-q+c,2)*dR(ki,2);

                % Compute recursively the entries of the
                % Hessian matrix
                % d^2 x1/du^2
                Hxxi(1,1) = Hxxi(1,1) + CP(i-p+b,j-q+c,1)*ddR(ki,1);
                % d^2 x1/dv^2
                Hxxi(1,2) = Hxxi(1,2) + CP(i-p+b,j-q+c,1)*ddR(ki,2);
                % d^2 x1/dudv
                Hxxi(1,3) = Hxxi(1,3) + CP(i-p+b,j-q+c,1)*ddR(ki,3);
                % d^2 x2/du^2
                Hxxi(2,1) = Hxxi(2,1) + CP(i-p+b,j-q+c,2)*ddR(ki,1);
                % d^2 x2/dv^2
                Hxxi(2,2) = Hxxi(2,2) + CP(i-p+b,j-q+c,2)*ddR(ki,2);
                % d^2 x2/dudv
                Hxxi(2,3) = Hxxi(2,3) + CP(i-p+b,j-q+c,2)*ddR(ki,3);
            end
        end

        % On the computation of the second derivatives of the 
        % basis functions in the physical space given those at
        % the parent domain

        % Initialize tranformation matrix for the second
        % derivatives
        DJxxi = zeros(5,5);

        % Assign the values of DJxxi
        % First block of DJxxi
        for ki=1:2
            for kj=1:2
                DJxxi(ki,kj) = Jxxi(ki,kj)^2;
            end
        end

        % Complete first row
        DJxxi(1,3) = 2*Jxxi(1,1)*Jxxi(1,2);
        DJxxi(1,4) = Hxxi(1,1)^2;
        DJxxi(1,5) = Hxxi(2,1)^2;

        % Complete second row
        DJxxi(2,3) = 2*Jxxi(2,1)*Jxxi(2,2);
        DJxxi(2,4) = Hxxi(1,2)^2;
        DJxxi(2,5) = Hxxi(2,2)^2;

        % Complete third row
        DJxxi(3,1) = Jxxi(1,1)*Jxxi(2,1);
        DJxxi(3,2) = Jxxi(1,2)*Jxxi(2,2);
        DJxxi(3,3) = Jxxi(1,1)*Jxxi(2,2)+Jxxi(2,1)*Jxxi(1,2);
        DJxxi(3,4) = Hxxi(1,3);
        DJxxi(3,5) = Hxxi(2,3);

        % Last block is the Jacobian matrix
        for ki=1:2
            for kj=1:2
                DJxxi(ki+3,kj+3) = Jxxi(ki,kj);
            end
        end

        % Compute the inverse of the Jacobian
        %invJxxi = inv(Jxxi);

        % The second order derivatives of the basis functions
        % with respect to the physical space

        % Initialization of the matrix with the derivatives
        dRdx = zeros(nNode_e,3);

        % loop over all the elements
        for ki=1:nNode_e
            deriv_R_actual = [ddR(ki,1:3) dR(ki,1:2)];
            deriv_Rdx = DJxxi\deriv_R_actual';
            dRdx(ki,1:3) = deriv_Rdx(1:3);
        end   

        % The element control variables vector affecting the 
        % current knot span
        d_actual(:,1) = d_el(i,j,:);

        % In the discrete form the residual equation yields
        % e = Ld*R*d_actual = Lp*C*L*Rs*d_actual, at first the
        % product Ld*R will be computed
        LdR = zeros(2,ndof_e);          

        % Loop over all the local degrees of freedom
        for ki=1:nNode_e
            % Assign the values of LdR recursively
            LdR(1,2*ki-1) = dRdx(ki,1)+dRdx(ki,2)*(1-parameters.nu)/2;
            LdR(1,2*ki) = dRdx(ki,3)*(1+parameters.nu)/2;
            LdR(2,2*ki-1) = dRdx(ki,3)*(1+parameters.nu)/2;
            LdR(2,2*ki) = dRdx(ki,2)+dRdx(ki,1)*(1-parameters.nu)/2;
        end

        % Multiply the result via the factor Eh/(1-nu^2) to get
        % the physical operator values
        LdR(:,:) = LdR(:,:)*parameters.E*parameters.t/(1-parameters.nu^2);

        % Compute the residual force on the Gauss point
        e_GP(k,l,1:2) = LdR*d_actual;
        
        k=k+1;
        u=u+s;
    end
    l=l+1;
    v=v+r;
end

figure(graph.index)

% First window, the x-component of the residual forces
subplot(2,1,1);

% Create the colored surface
surf(P(:,:,1),P(:,:,2),P(:,:,3),e_GP(:,:,1));

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
title ('x-component of the residual forces');

% Second window, the x-component of the residual forces
subplot(2,1,2);

% Create the colored surface
surf(P(:,:,1),P(:,:,2),P(:,:,3),e_GP(:,:,2));

% Graphics options
shading interp;
colormap('default');

% invert default colormap => red=negativ, blue=positive
COL=colormap;
invCOL(:,1)=COL(:,3);
invCOL(:,2)=COL(:,2);
invCOL(:,3)=COL(:,1);
colormap(invCOL);

% make colormap symmetric
colim = caxis;
caxis([-max(abs(colim)) max(abs(colim))]);

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
title ('y-component of the residual forces');

hold off;

% Update the graph index
index = graph.index + 1;

end

