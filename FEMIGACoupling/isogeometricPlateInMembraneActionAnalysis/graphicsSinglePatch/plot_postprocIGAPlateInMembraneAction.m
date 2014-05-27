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
function index = plot_postprocIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,rb,parameters,Fl,dHat,graph,outMsg)
%% Function documentation
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

% Initialize handle to the figure
figure(graph.index)

% Grid point number for the plotting of both the B-Spline surface the knots
% as well as the resultant computation over the domain
xiGrid = 49;
etaGrid = 49;

%% 1. Plot the first window: Reference and/or current configuration

% Plot the window
subplot(2,1,1);
plot_postprocCurrentConfigurationIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,xiGrid,etaGrid,rb,Fl,dHat,graph);

% Assign graphic properties and title
camlight left; lighting phong;
view(2);
axis equal;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
if strcmp(graph.postprocConfig,'reference')
    title('Reference configuration');
elseif strcmp(graph.postprocConfig,'current')
    title('Current configuration');
elseif strcmp(graph.postprocConfig,'referenceCurrent')
    title('Reference and current configuration');
end
hold off;

%% 2. Plot the second window: Resultant visualization

% Plot the window
subplot(2,1,2);
plot_postprocResultantsIGAPlateInMembraneAction(p,q,Xi,Eta,CP,isNURBS,parameters,xiGrid,etaGrid,dHat,graph);

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
colorbar;

view(2);
axis equal;
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