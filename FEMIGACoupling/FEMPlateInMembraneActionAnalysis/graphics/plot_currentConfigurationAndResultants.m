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
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = plot_currentConfigurationAndResultants(mesh,rb,displacement,parameters,analysis,resultant,component,graph)
%% Function documentation
%
% Plots the current configuration of a plate in membrane action given the
% displacement field and visualizes the displacement/strain or stress field
% over the initial configuration of the plate
%
%              input :
%               mesh : Elements and nodes of the mesh
%                 rb : Vector of the Dirichlet boundary conditions with their
%                      global numbering
%       displacement : The displacement field sorted in a vector according to its
%                      global numbering
%         parameters : The material properties of the structure
%           analysis : Analysis type (plane stress or plane strain)
%          resultant : Resultant to be visualized
%          component : component of the resultant to be visualized
%              graph : Structure on the graphics (indices etc.)
%
%             output :
%              index : The index of the current graph
%
%% Function main body

%% 0. Read input

% Number of DoFs at the element level (depends on the element type)
noNodesElement = 3;
noDoFsElement = noNodesElement*2;

%% 1. Compute the new loactions for the vertices of the triangles in the mesh

% Initialize the array of the displaced nodes
nodes_displaced = zeros(length(mesh.nodes),3);

% Initialize pseudocounter
counter = 1;

for i=1:length(mesh.nodes)
    % Add the x and y components of the displacement field
    nodes_displaced(i,1) = mesh.nodes(i,1) + displacement(2*counter-1);
    nodes_displaced(i,2) = mesh.nodes(i,2) + displacement(2*counter);
    
    % Update counter
    counter = counter + 1;
end

% Initialize figure
figure(graph.index)

%% 2. 1st window: the deformed and the undeformed configuration
subplot(2,1,1);

% Visualize the displaced elements on the mesh

% Reference configuration
if strcmp(graph.visualization.geometry,'reference_and_current')||strcmp(graph.visualization.geometry,'current');
    patch('faces',mesh.elements,'vertices',nodes_displaced(:,1:2),'facecolor','g','edgecolor','black');
    hold on;
end
% Current configuration
if strcmp(graph.visualization.geometry,'reference_and_current')||strcmp(graph.visualization.geometry,'reference');
    patch('faces',mesh.elements,'vertices',mesh.nodes(:,1:2),'facecolor','none','edgecolor','black');
end
axis equal off;
axis on;
title('The current configuration of the mesh');

% Visualize the Dirichlet boundary conditions on the mesh

% Create the supports
[xs,ys,zs] = createSupports(nodes_displaced,rb);

%supports
hold on;
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end
hold off;

% Graph properties
view (2);
axis equal;
axis on;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);

% Title
if strcmp(analysis.type,'PLANE_STESS')
    title('Deformed configuration corresponding to plain stress analysis');
elseif strcmp(analysis.type,'PLANE_STRAIN')
    title('Deformed configuration corresponding to plain strain analysis');
end

%% 3. 2nd window: Visualization of the resultants
subplot(2,1,2);

% Grid at each triangular element
grid_xi = 5;
grid_eta = 5;

% step size to both directions
step_xi = 1/(grid_xi+1);
step_eta = 1/(grid_eta+1);

% Plot the mesh edges
plot_EdgesTriangularMesh2D(mesh);
hold on;

for c_element=1:size(mesh.elements(:,1))
    % Compute the basis functions evaluated at the grid points
    P = zeros(grid_xi+2,grid_eta+2,3);
    tensorialQuantityOverDomain = zeros(grid_xi+2,grid_eta+2);

    % Initialize local coordinates
    xi = 0;
    eta = 0;
    
    % get element from the mesh
    element = mesh.elements(c_element,:);
    
    % get coordinates of the element vertices
    nodes = mesh.nodes(element,:);
    node_i = element(1,1);
    node_j = element(1,2);
    node_k = element(1,3);
    
    % The vertices of the current triangle
    Pi = mesh.nodes(node_i,:);
    Pj = mesh.nodes(node_j,:);
    Pk = mesh.nodes(node_k,:);
    
    % Element freedom table
    EFT = zeros(noDoFsElement,1);
    for j=1:noNodesElement
        EFT(2*j-1,1) = 2*element(j)-1;
        EFT(2*j,1) = 2*element(j);
    end
    
    % Get the element displacement vector
    displacementElement = displacement(EFT);

    % The moving vertices
    Pi_eta = Pi;
    Pj_eta = Pj;
    
    % Loop over all the sampling points
    for j=1:grid_eta+2
        for i=1:grid_xi+2
            
            % Initializations for the velocities
            ux = 0;
            uy = 0;
            
            % Get the point in the interior of the line defined from Pi and Pj
            P(i,j,:) = xi*Pi_eta+(1-xi)*Pj_eta;
        
            % Evaluate the CST basis functions at P
            if strcmp(resultant,'displacement')
                N = computeCST2DBasisFunctions(Pi,Pj,Pk,P(i,j,1),P(i,j,2));
            elseif strcmp(resultant,'strain') || strcmp(resultant,'stress')
                [dN,~] = computeCST2DBasisFunctionsAndFirstDerivatives(Pi,Pj,Pk,P(i,j,1),P(i,j,2));
            end
            
            % Evaluate the tensorial field value on P by looping over all 
            % the products of the basis functions with the nodal solution
            
            if strcmp(resultant,'displacement')
                % Compute the displacement field
                for k = 1:noNodesElement
                    % Displacement in x-direction
                    ux = ux + N(k) * displacementElement(2*k-1);
                
                    % Displacement in y-direction
                    uy = uy + N(k) * displacementElement(2*k);
                end
                if strcmp(component,'x')
                    tensorialQuantityOverDomain(i,j) = ux;
                elseif strcmp(component,'y')
                    tensorialQuantityOverDomain(i,j) = uy;
                elseif strcmp(component,'2norm')
                    tensorialQuantityOverDomain(i,j) = sqrt(ux^2+uy^2);
                end
            elseif strcmp(resultant,'strain')
                % Compute the B-Operator Matrix    
                B = [dN(1,2) 0       dN(2,2) 0       dN(3,2) 0
                     0       dN(1,3) 0       dN(2,3) 0       dN(3,3)
                     dN(1,3) dN(1,2) dN(2,3) dN(2,2) dN(3,3) dN(3,2)];
                    
                % Compute the strain field in Voigt notation
                strain = B*displacementElement;
                if strcmp(component,'x')
                    tensorialQuantityOverDomain(i,j) = strain(1);
                elseif strcmp(component,'y')
                    tensorialQuantityOverDomain(i,j) = strain(2);
                elseif strcmp(component,'xy')
                    tensorialQuantityOverDomain(i,j) = strain(3);
                end
            elseif strcmp(resultant,'stress')
                % Compute the B-Operator Matrix    
                B = [dN(1,2) 0       dN(2,2) 0       dN(3,2) 0
                     0       dN(1,3) 0       dN(2,3) 0       dN(3,3)
                     dN(1,3) dN(1,2) dN(2,3) dN(2,2) dN(3,3) dN(3,2)];
                    
                % Compute the strain field in Voigt notation
                strain = B*displacementElement;
                
                % Compute the material matrix
                if strcmp(analysis.type,'PLANE_STRESS')
                    preFactor = parameters.E/(1-parameters.nue^2);
                    C = preFactor*[1                parameters.nue 0
                                   parameters.nue   1              0
                                   0                0              (1-parameters.nue)/2];
                elseif strcmp(analysis.type,'PLAIN_STRAIN')
                    preFactor = parameters.E*(1-parameters.nue)/(1+parameters.nue)/(1-2*parameters.nue);
                    C = preFactor*[1                                 parameters.nue/(1-parameters.nue) 0
                                   parameters.nue/(1-parameters.nue) 1                                 0
                                   0                                 0                                (1-2*parameters.nue)/2/(1-parameters.nue)];
                end
                
                % compute the stress field in Voigt notation
                stress = C*strain;
                if strcmp(component,'x')
                    tensorialQuantityOverDomain(i,j) = stress(1);
                elseif strcmp(component,'y')
                    tensorialQuantityOverDomain(i,j) = stress(2);
                elseif strcmp(component,'xy')
                    tensorialQuantityOverDomain(i,j) = stress(3);
                end
            end
                  
            % Update local coordinate by the step size
            xi = xi + step_xi;
        end
        % Reset xi local coordinate
        xi = 0;
    
    	% Update eta local coordinate
        eta = eta + step_eta;
    
        % Update the moving vertices over the triangle edges
        Pi_eta = eta*Pk+(1-eta)*Pi;
        Pj_eta = eta*Pk+(1-eta)*Pj;
    end


    % Graphical output  of the basis functions
    surf(P(:,:,1),P(:,:,2),P(:,:,3),tensorialQuantityOverDomain(:,:,1),'EdgeColor','none');
    hold on;
end

% Graphics options
colormap('default');
colorbar;
hold on;

% graph properties
view (2);
axis equal;
axis on;
grid on;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if strcmp(resultant,'displacement')
    if strcmp(component,'x')
        title('x-component of the displacement field u_x');
    elseif strcmp(component,'y')
        title('y-component of the displacement field u_y');
    elseif strcmp(component,'2norm')
        title('2-norm of the displacement field ||u||_2');
    end
elseif strcmp(resultant,'strain')
    if strcmp(component,'x')
        title('normal xx-component of the strain field epsilon_{xx}');
    elseif strcmp(component,'y')
        title('normal yy-component of the strain field epsilon_{yy}');
    elseif strcmp(component,'xy')
        title('shear xy-component of the strain field epsilon_{xy}');
    end
elseif strcmp(resultant,'stress')
    if strcmp(component,'x')
        title('normal xx-component of the stress field sigma_{xx}');
    elseif strcmp(component,'y')
        title('normal yy-component of the stress field sigma_{yy}');
    elseif strcmp(component,'xy')
        title('shear xy-component of the stress field sigma_{xy}');
    end
end

hold off;

% Update index
index = graph.index + 1;

end

