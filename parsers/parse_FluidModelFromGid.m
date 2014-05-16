%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   Aditya Ghantasala                  (aditya.ghantasala@tum.de)         %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fldMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,analysis,parameters,nLinearAnalysis,fldDynamics,gaussInt] = ...
    parse_FluidModelFromGid(pathToCase,caseName,outMsg)
%% Function documentation
%
% Parses data from an input file created using GiD for a fluid boundary 
% value problem.
%
%           Input :
%      pathToCase : The absolute path to the inputGiD case folder
%        caseName : The name of the case in the inputGiD case folder
%          outMsg : On the output information on the command window
%
%          Output :
%          fldMsh :     .nodes : The nodes in the FE mesh
%                    .elements : The elements in the FE mesh
%          homDBC : The global numbering of the nodes where homogeneous
%                   Dirichlet boundary conditions are applied
%        inhomDBC : The global numbering of the nodes where inhomogeneous
%                   Dirichlet boundary conditions are applied
%  valuesInhomDBC : The prescribed values for the inhomogeneous Dirichlet
%                   boundary conditions
%             NBC :     .nodes : The nodes where Neumann boundary 
%                                conditions are applied
%                    .loadType : The type of the load for each Neumann node
%                   .fctHandle : The function handle for each Neumann node
%                                for the computation of the load vector 
%                                (these functions are unde the folder load)
%        analysis : .type : The analysis type
%      parameters : Problem specific technical parameters
% nLinearAnalysis :     .scheme : The employed nonlinear scheme
%                    .tolerance : The residual tolerance
%                      .maxIter : The maximum number of the nonlinear 
%                                 iterations
%     fldDynamics : .timeDependence : Steady-state or transient analysis
%                           .scheme : The time integration scheme
%                               .T0 : The start time of the simulation
%                             .TEnd : The end time of the simulation
%                          .nTSteps : The number of the time steps
%        gaussInt : On the Gauss Point integration
%                             .type : 'default', 'user'
%                       .domainNoGP : Number of Gauss Points for the domain
%                                     integration
%                     .boundaryNoGP : Number of Gauss Points for the
%                                     boundary integration
% Function layout :
%
% 1. Load the input file from GiD
%
% 2. Load the analysis type
%
% 3. Load the material properties
%
% 4. Load the nonlinear scheme
%
% 5. Load the time integration scheme
%
% 6. Load the Gauss Point integration scheme
%
% 7. Load the structural nodes
%
% 8. Load the structural elements by connectivity arrays
%
% 9. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
%
% 10. Load the nodes on the Neumann boundary together with the load application information
%
% 11. Get edge connectivity arrays for the Neumann edges
%
% 12. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________\n');
    fprintf('###########################################################\n');
    fprintf('Parsing data from GiD input file for a fluid boundary value\n');
    fprintf('problem has been initiated\n');
    fprintf('___________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Initialize output arrays
homDBC = [];
inhomDBC = [];
valuesInhomDBC = [];

%% 1. Load the input file from GiD
fstring = fileread([pathToCase caseName '.dat']); 

%% 2. Load the analysis type
block = regexp(fstring,'FLUID_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
analysis.type = out{1}{2};
fprintf('>> Analysis type: %s \n',analysis.type);

%% 3. Load the material properties
block = regexp(fstring,'FLUID_MATERIAL_PROPERTIES','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
parameters.rho = str2double(out{1}{2});
parameters.nue = str2double(out{1}{4});

%% 4. Load the nonlinear scheme
block = regexp(fstring,'FLUID_NLINEAR_SCHEME','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
nLinearAnalysis.scheme = out{1}{2};
nLinearAnalysis.tolerance = str2double(out{1}{4});
nLinearAnalysis.maxIter = str2double(out{1}{6});

%% 5. Load the time integration scheme
block = regexp(fstring,'FLUID_TRANSIENT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',' ','MultipleDelimsAsOne', 1);
fldDynamics.timeDependence = out{1}{2};
fldDynamics.scheme = out{1}{4};
if strcmp(fldDynamics.scheme,'bossak')
    fldDynamics.alphaBeta =  str2double(out{1}{6});
    fldDynamics.gamma =  str2double(out{1}{8});
end
fldDynamics.T0 = str2double(out{1}{10});
fldDynamics.TEnd = str2double(out{1}{12});
fldDynamics.noTimeSteps = str2double(out{1}{14});
fldDynamics.isAdaptive = out{1}{16};
fprintf('>> Fluid dynamics: %s \n',fldDynamics.timeDependence);
if ~strcmp(fldDynamics.timeDependence,'STEADY_STATE')
    fprintf('\t>> Time integration scheme: %s \n',fldDynamics.scheme);
    if strcmp(fldDynamics.scheme,'bossak')
        fprintf('\t \t>> alphaBeta =  %s \n',fldDynamics.alphaBeta);
        fprintf('\t \t>> gamma =  %s \n',fldDynamics.gamma);
    end
    fprintf('\t>> Start time of the simulation: %f \n',fldDynamics.T0);
    fprintf('\t>> End time of the simulation: %f \n',fldDynamics.TEnd);
    fprintf('\t>> Number of time steps: %f \n',fldDynamics.noTimeSteps);
    fldDynamics.dt = (fldDynamics.TEnd - fldDynamics.T0)/fldDynamics.noTimeSteps;
    fprintf('\t>> Time step size: %f \n',fldDynamics.dt);
end

%% 6. Load the Gauss Point integration scheme
block = regexp(fstring,'FLUID_INTEGRATION','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',' ','MultipleDelimsAsOne', 1);
gaussInt.type = out{1}{2};
if strcmp(gaussInt.type,'user')
    gaussInt.domainNoGP = str2double(out{1}{4});
    gaussInt.boundaryNoGP = str2double(out{1}{6});
end
fprintf('>> Gauss integration type: %s \n',gaussInt.type);
if strcmp(gaussInt.type,'USER')
    fprintf('\t>> No. Gauss Points for the domain integration: %d \n',gaussInt.domainNoGP);
    fprintf('\t>> No. Gauss Points for the boundary integration: %d \n',gaussInt.boundaryNoGP);
end

%% 7. Load the fluid nodes
block = regexp(fstring,'FLUID_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
fldMsh.nodes = out(:,2:4);
fprintf('>> Number of nodes in the mesh: %d \n',length(fldMsh.nodes(:,1)));

%% 8. Load the fluid elements by connectivity arrays
block = regexp(fstring,'FLUID_ELEMENTS','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
fldMsh.elements = out(:,2:4);
fprintf('>> Number of elements in the mesh: %d \n',length(fldMsh.elements));

%% 9. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
block = regexp(fstring,'FLUID_DIRICHLET_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);

% Filter out the actual DOFs
nDOFsPerNodeFromGiD = 4;
if strcmp(analysis.type,'NAVIER_STOKES')
    nDOFsPerNode = 3;
    nDOFsPerNodeArray = [1 2 4];
end

% Initialize counter
counterHomDBC = 1;
counterInhomDBC = 1;

% Find the number of nodes where Dirichlet boundary conditions are applied
nDBCNodes = length(out)/(nDOFsPerNodeFromGiD+1);

for i=1:nDBCNodes
    % Get the Dirichlet node ID
    nodeID = out((nDOFsPerNodeFromGiD+1)*i-nDOFsPerNodeFromGiD);
    
    % Get the x-component of the prescribed value
    for jCounter = 1:length(nDOFsPerNodeArray)
        % Get the actual j-counter
        j = nDOFsPerNodeArray(jCounter);
        
        % Get the prescribed value at the current DOF
        presValue = out((nDOFsPerNodeFromGiD+1)*i-nDOFsPerNodeFromGiD+j);
        if ~isnan(presValue)
            if presValue == 0
                homDBC(counterHomDBC) = nDOFsPerNode*nodeID-nDOFsPerNode+j;
                counterHomDBC = counterHomDBC + 1;
            else
                inhomDBC(counterInhomDBC) = nDOFsPerNode*nodeID-nDOFsPerNode+j;
                valuesInhomDBC(counterInhomDBC) = presValue;
                counterInhomDBC = counterInhomDBC + 1;
            end
        end
    end
end

% Sort out the vectors
homDBC = sort(homDBC);
[inhomDBC,indexSorting] = sort(inhomDBC);
valuesInhomDBC = valuesInhomDBC(indexSorting);

%% 10. Load the nodes on the Neumann boundary together with the load application information
block = regexp(fstring,'FLUID_FORCE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %s %s','delimiter',' ','MultipleDelimsAsOne', 1);
end

if ~isempty(out)
    out = out{1};
    NBC.nodes = cell2mat(out(:,1));
    outLoadType = out(:,2);
    NBC.loadType = cell2mat(outLoadType{1});
    outFctHandle = out(:,3);
    NBC.fctHandle = cell2mat(outFctHandle{1});
end

%% 11. Get edge connectivity arrays for the Neumann edges
if ~isempty(out)
    fprintf('>> Neumann boundary edges: %d \n',length(NBC.nodes)-1);

    % Initialize the Neumann boundary lines
    NBC.lines = zeros(length(NBC.nodes)-1,3);

    % Initialize line counter
    counterLines = 1;

    % Loop over each node pair
    for i=1:length(NBC.nodes)
        for j=i:length(NBC.nodes)
            % If we are not in the same node
            if i~=j
                % Get the node index in the element array
                nodeI = NBC.nodes(i);
                nodeJ = NBC.nodes(j);


                % Find the element indices to which the nodes belong
                [indexI,~] = find(nodeI == fldMsh.elements);
                [indexJ,~] = find(nodeJ == fldMsh.elements);

                % For all the element indices to which indexJ belongs to
                for k=1:length(indexJ)
                    % Find the common elements to which both nodes belong to
                    commonElmnts = find(indexJ(k) == indexI);

                    % If there are commont elements to which the nodes belong
                    if norm(commonElmnts)~=0
                        % Get the common element index
                        elementIndex = indexI(commonElmnts);

                        % Store the line into the NBC.line array with the same
                        % ordering as the are stored in the element array
                        NBC.lines(counterLines,:) = ...
                            [NBC.nodes(i) NBC.nodes(j) elementIndex];

                        % Update counter
                        counterLines = counterLines + 1;
                    end
                end
            end
        end
    end
else
    NBC = 'undefined';
end

%% 12. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nParsing took %.2d seconds \n\n',computationalTime);
    fprintf('_______________________Parsing Ended_______________________\n');
    fprintf('###########################################################\n\n\n');
end

end