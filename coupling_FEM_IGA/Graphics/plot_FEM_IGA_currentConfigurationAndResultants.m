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
function index = plot_FEM_IGA_currentConfigurationAndResultants(p,q,Xi,Eta,CP,isNURBS,rb,parameters,Fl,dHatIGA,graph,outMsg,...
    mesh,rbFEM,dHatFEM,parametersFEM,analysis)
%% Function documentation
%
% Plots the current configuration of a plate in membrane action given the
% displacement field and visualizes the displacement/strain or stress field
% over the initial configuration of the plate
%
%
%%%%%%%%%%%%%           IGA part
%
% Plots two windows in one: The first contains the reference and/or the
% current configuration of an isogeometric plate in membrane action given
% the Control Point displacement. The second window contains the
% visualization of the selected resultant component over the plate's domain
% in the reference configuration.
%
%       Input :
%         p,q : Polynomial degrees
%      Xi,Eta : Knot vectors in xi,eta-direction
%          CP : Control point coordinates and weights of the undeformed 
%               plate
%     isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
%          rb : Vector containing information on the supports
%  parameters : Technical and geometrical parameters of the plate
%          Fl : The applied load vector
%        dHat : The displacement field of the control points
%       graph : Information on the graphics
%      outMsg : Whether or not to output message on refinement progress
%               'outputEnabled' : enables output information
%
%      Output :
%       index : The index of the current graph
%               
%               graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Plot the first window: Reference and/or current configuration
%
% 2. Plot the second window: Resultant visualization
%
% 3. Update the graph index
%
% 4. Appendix
%
%
%
%%%%%%%%%%%%           FEM part
%
%
%              input :
%               mesh : Elements and nodes of the mesh
%              rbFEM : Vector of the Dirichlet boundary conditions with their
%                      global numbering
%       displacement : The displacement field sorted in a vector according to its
%                      global numbering
%      parametersFEM : The material properties of the structure
%           analysis : Analysis type (plane stress or plane strain)
%         
%         
%
%
%% Function main body


%%%%%%%%%%%%%%%%%% IGA part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function documentation
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________________\n');
    fprintf('##################################################################\n');
    fprintf('Plotting postprocessing configuration and resultant visualization\n');
    fprintf('for the isogeometric plate in mebrane action has been initiated\n\n');
    fprintf('Configuration to be visualized (1st window): ');
    if strcmp(graph.postprocConfig,'reference')
        fprintf('Reference\n');
    elseif strcmp(graph.postprocConfig,'current')
        fprintf('Current\n');
    elseif strcmp(graph.postprocConfig,'referenceCurrent')
        fprintf('Reference and current\n');
    end
    fprintf('Resultant to be visualized (2nd window): ');
    if strcmp(graph.resultant,'displacement')
        if strcmp(graph.component,'x')
            fprintf('Displacement component u_x\n');
        elseif strcmp(graph.component,'y')
            fprintf('Displacement component u_y\n');
        elseif strcmp(graph.component,'xy')
            fprintf('Displacement magnitude ||u||_2\n');
        end
    elseif strcmp(graph.resultant,'strain')
        if strcmp(graph.component,'x')
            fprintf('Strain component epsilon_xx\n');
        elseif strcmp(graph.component,'y')
            fprintf('Strain component epsilon_yy\n');
        elseif strcmp(graph.component,'xy')
            fprintf('Strain component epsilon_xy\n');
        end
    elseif strcmp(graph.resultant,'stress')
        if strcmp(graph.component,'x')
            fprintf('Stress component sigma_xx\n');
        elseif strcmp(graph.component,'y')
            fprintf('Stress component sigma_yy\n');
        elseif strcmp(graph.component,'xy')
            fprintf('Stress component sigma_xy\n');
        end
    end
    fprintf('__________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input


%%%%%%%%%%%%%%%%%%%  IGA

% Grid point number for the plotting of both the B-Spline surface the knots
% as well as the resultant computation over the domain
xiGrid = 49;
etaGrid = 49;


%%%%%%%%%%%%%%%%%%  FEM


% Number of DoFs at the element level (depends on the element type)
noNodesElement = 3;
noDoFsElement = noNodesElement*2;

%%%% Compute the new loactions for the vertices of the triangles in the mesh

% Initialize the array of the displaced nodes
nodes_displaced = zeros(length(mesh.nodes),3);

% Initialize pseudocounter
counter = 1;

for i=1:length(mesh.nodes)
    % Add the x and y components of the displacement field
    nodes_displaced(i,1) = mesh.nodes(i,1) + dHatFEM(2*counter-1);
    nodes_displaced(i,2) = mesh.nodes(i,2) + dHatFEM(2*counter);
    
    % Update counter
    counter = counter + 1;
end

%%%%%%%%%%%%%%%%%%  General

% Initialize figure
figure(graph.index)



%% 1. Plot the first window: Reference and/or current configuration

%%%%%%%%%%%%%%%%%%%  IGA


% Plot the window
subplot(2,1,1);
plot_postprocCurrentConfigurationIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,xiGrid,etaGrid,rb,Fl,dHatIGA,graph);

% Assign graphic properties and title
camlight left; lighting phong;
view(2);
if strcmp(graph.postprocConfig,'reference')
    title('Reference configuration');
elseif strcmp(graph.postprocConfig,'current')
    title('Current configuration');
elseif strcmp(graph.postprocConfig,'referenceCurrent')
    title('Reference and current configuration');
end
hold off;

%%%%%%%%%%%%%%%%%%  FEM

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

% Visualize the Dirichlet boundary conditions on the mesh

% Create the supports
[xs,ys,zs] = createSupports(nodes_displaced,rbFEM);

%supports
for k =1:length(xs(:,1))
    plot3(xs(k,:),ys(k,:),zs(k,:),'Linewidth',2,'Color','black');
end
hold off;

% Graph properties
view (2);
% axis ([0 6 0 6])

% Title
if strcmp(analysis.type,'PLANE_STESS')
    title('Deformed configuration corresponding to plain stress analysis');
elseif strcmp(analysis.type,'PLANE_STRAIN')
    title('Deformed configuration corresponding to plain strain analysis');
end

axis tight;
axis fill;
axis on;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);

%% 2. Plot the second window: Resultant visualization

%%%%%%%%%%%%%%%%%%%  IGA


% Plot the window
subplot(2,1,2);
plot_postprocResultantsIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,parameters,xiGrid,etaGrid,dHatIGA,graph);

% Assign graphic properties and title
if strcmp(graph.resultant,'displacement')
    if strcmp(graph.component,'x')
        titleString = 'Displacement component d_x';
    elseif strcmp(graph.component,'y')
        titleString = 'Displacement component d_y';
    elseif strcmp(graph.component,'2norm')
        titleString = 'Displacement magnitude ||d||_2';
    end
elseif strcmp(graph.resultant,'strain')
    if strcmp(graph.component,'x')
        titleString = 'Strain component \epsilon_{xx}';
    elseif strcmp(graph.component,'y')
        titleString = 'Strain component \epsilon_{yy}';
    elseif strcmp(graph.component,'xy')
        titleString = 'Strain component \epsilon_{xy}';
    elseif strcmp(graph.component,'1Principal')
        titleString = '1st pricipal strain field \epsilon_1';
    elseif strcmp(graph.component,'2Principal')
        titleString = '2nd pricipal strain field \epsilon_2';
    end
elseif strcmp(graph.resultant,'stress')
    if strcmp(graph.component,'x')
        titleString = 'Stress component \sigma_{xx}';
    elseif strcmp(graph.component,'y')
        titleString = 'Stress component \sigma_{yy}';
    elseif strcmp(graph.component,'xy')
        titleString = 'Stress component \sigma_{xy}';
    elseif strcmp(graph.component,'1Principal')
        titleString = '1st pricipal stress field \sigma_1';
    elseif strcmp(graph.component,'2Principal')
        titleString = '2nd pricipal stress field \sigma_2';
    end
end
title(titleString);

shading interp;
colormap('default');

% invert default colormap => red = negativ, blue = positive
% COL = colormap;
% invCOL(:,1) = COL(:,3);
% invCOL(:,2) = COL(:,2);
% invCOL(:,3) = COL(:,1);
% colormap(invCOL);

% make colormap symmetric
% colim = caxis;
% caxis([-max(abs(colim)) max(abs(colim))]);



%%%%%%%%%%%%%%%%%%  FEM


% Grid at each triangular element
grid_xi = 5;
grid_eta = 5;

% step size to both directions
step_xi = 1/(grid_xi+1);
step_eta = 1/(grid_eta+1);

% Plot the mesh edges
hold on;
plot_EdgesTriangularMesh2D(mesh);


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
    displacementElement = dHatFEM(EFT);

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
            if strcmp(graph.resultant,'displacement')
                N = computeCST2DBasisFunctions(Pi,Pj,Pk,P(i,j,1),P(i,j,2));
            elseif strcmp(graph.resultant,'strain') || strcmp(graph.resultant,'stress')
                [dN,~] = computeCST2DBasisFunctionsAndFirstDerivatives(Pi,Pj,Pk,P(i,j,1),P(i,j,2));
            end
            
            % Evaluate the tensorial field value on P by looping over all 
            % the products of the basis functions with the nodal solution
            
            if strcmp(graph.resultant,'displacement')
                % Compute the displacement field
                for k = 1:noNodesElement
                    % Displacement in x-direction
                    ux = ux + N(k) * displacementElement(2*k-1);
                
                    % Displacement in y-direction
                    uy = uy + N(k) * displacementElement(2*k);
                end
                if strcmp(graph.component,'x')
                    tensorialQuantityOverDomain(i,j) = ux;
                elseif strcmp(graph.component,'y')
                    tensorialQuantityOverDomain(i,j) = uy;
                elseif strcmp(graph.component,'2norm')
                    tensorialQuantityOverDomain(i,j) = sqrt(ux^2+uy^2);
                end
            elseif strcmp(graph.resultant,'strain')
                % Compute the B-Operator Matrix    
                B = [dN(1,2) 0       dN(2,2) 0       dN(3,2) 0
                     0       dN(1,3) 0       dN(2,3) 0       dN(3,3)
                     dN(1,3) dN(1,2) dN(2,3) dN(2,2) dN(3,3) dN(3,2)];
                    
                % Compute the strain field in Voigt notation
                strain = B*displacementElement;
                if strcmp(graph.component,'x')
                    tensorialQuantityOverDomain(i,j) = strain(1);
                elseif strcmp(graph.component,'y')
                    tensorialQuantityOverDomain(i,j) = strain(2);
                elseif strcmp(graph.component,'xy')
                    tensorialQuantityOverDomain(i,j) = strain(3);
                end
            elseif strcmp(graph.resultant,'stress')
                % Compute the B-Operator Matrix    
                B = [dN(1,2) 0       dN(2,2) 0       dN(3,2) 0
                     0       dN(1,3) 0       dN(2,3) 0       dN(3,3)
                     dN(1,3) dN(1,2) dN(2,3) dN(2,2) dN(3,3) dN(3,2)];
                    
                % Compute the strain field in Voigt notation
                strain = B*displacementElement;
                
                % Compute the material matrix
                if strcmp(analysis.type,'PLANE_STRESS')
                    preFactor = parametersFEM.E/(1-parametersFEM.nue^2);
                    C = preFactor*[1                parametersFEM.nue 0
                                   parametersFEM.nue   1              0
                                   0                0              (1-parametersFEM.nue)/2];
                elseif strcmp(analysis.type,'PLAIN_STRAIN')
                    preFactor = parametersFEM.E*(1-parametersFEM.nue)/(1+parametersFEM.nue)/(1-2*parametersFEM.nue);
                    C = preFactor*[1                                 parametersFEM.nue/(1-parametersFEM.nue) 0
                                   parametersFEM.nue/(1-parametersFEM.nue) 1                                 0
                                   0                                 0                                (1-2*parametersFEM.nue)/2/(1-parametersFEM.nue)];
                end
                
                % compute the stress field in Voigt notation
                stress = C*strain;
                if strcmp(graph.component,'x')
                    tensorialQuantityOverDomain(i,j) = stress(1);
                elseif strcmp(graph.component,'y')
                    tensorialQuantityOverDomain(i,j) = stress(2);
                elseif strcmp(graph.component,'xy')
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
% axis auto;
% axis ([0 6 0 6])

axis tight;
axis fill;
axis on;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);


%% 3. Update the graph index
index = graph.index + 1;

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Plotting the current configuration took %.2d seconds \n\n',computationalTime);
    fprintf('______________Plotting Current Configuration Ended________________\n');
    fprintf('##################################################################\n\n\n');
end

end