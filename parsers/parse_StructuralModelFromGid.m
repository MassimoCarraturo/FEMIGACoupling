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
function [strMsh,homDBC,inhomDBC,valuesInhomDBC,NBC,IBC,analysis,parameters,nLinearAnalysis,strDynamics] = ...
    parse_StructuralModelFromGid(pathToCase,caseName,outMsg)
%% Function documentation
%
% Parses data from an input file created using GiD for a structural 
% boundary value problem.
%
%           Input :
%      pathToCase : The absolute path to the inputGiD case folder
%        caseName : The name of the case in the inputGiD case folder
%          outMsg : On the output information on the command window
%
%          Output :
%          strMsh :     .nodes : The nodes in the FE mesh
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
%             IBC :     .nodes : The nodes on the interface boundary 
%        analysis : .type : The analysis type
%      parameters : Problem specific technical parameters
% nLinearAnalysis :     .scheme : The employed nonlinear scheme
%                    .tolerance : The residual tolerance
%                      .maxIter : The maximum number of the nonlinear 
%                                 iterations
%     strDynamics : .timeDependence : Steady-state or transient analysis
%                           .scheme : The time integration scheme
%                               .T0 : The start time of the simulation
%                             .TEnd : The end time of the simulation
%                          .nTSteps : The number of the time steps
%
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
% 6. Load the structural nodes
%
% 7. Load the structural elements by connectivity arrays
%
% 8. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
%
% 9. Load the nodes on the Neumann boundary together with the load application information
%
% 10. Get edge connectivity arrays for the Neumann edges
%
% 11. Load the nodes on the IGA/FEM interface
%
% 12. Get edge connectivity arrays for the IGA/FEM interface
%
% 13. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('________________________________________________________________\n');
    fprintf('################################################################\n');
    fprintf('Parsing data from GiD input file for a structural boundary value\n');
    fprintf('problem has been initiated\n');
    fprintf('________________________________________________________________\n\n');

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
block = regexp(fstring,'STRUCTURE_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
analysis.type = out{1}{2};
fprintf('>> Analysis type: %s \n',analysis.type);

%% 3. Load the material properties
block = regexp(fstring,'STRUCTURE_MATERIAL_PROPERTIES','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
parameters.rho = str2double(out{1}{2});
parameters.E = str2double(out{1}{4});
parameters.nue = str2double(out{1}{6});

%% 4. Load the nonlinear scheme
block = regexp(fstring,'STRUCTURE_NLINEAR_SCHEME','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',',','MultipleDelimsAsOne', 1);
nLinearAnalysis.scheme = out{1}{2};
nLinearAnalysis.tolerance = str2double(out{1}{4});
nLinearAnalysis.maxIter = str2double(out{1}{6});

%% 5. Load the time integration scheme
block = regexp(fstring,'STRUCTURE_TRANSIENT_ANALYSIS','split');
block(1) = [];
out = textscan(block{1},'%s','delimiter',' ','MultipleDelimsAsOne', 1);
strDynamics.timeDependence = out{1}{2};
strDynamics.scheme = out{1}{4};
strDynamics.T0 = str2double(out{1}{6});
strDynamics.TEnd = str2double(out{1}{8});
strDynamics.nTSteps = str2double(out{1}{10});
fprintf('>> Structural dynamics: %s \n',strDynamics.timeDependence);
if ~strcmp(strDynamics.timeDependence,'STEADY-STATE')
    fprintf('\t>> Time integration scheme: %s \n',strDynamics.scheme);
    fprintf('\t>> Start time of the simulation: %f \n',strDynamics.T0);
    fprintf('\t>> End time of the simulation: %f \n',strDynamics.TEnd);
    fprintf('\t>> Number of time steps: %f \n',strDynamics.nTSteps);
end

%% 6. Load the structural nodes
block = regexp(fstring,'STRUCTURE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
strMsh.nodes = out(:,2:4);
fprintf('>> Number of nodes in the mesh: %d \n',length(strMsh.nodes(:,1)));

%% 7. Load the structural elements by connectivity arrays
block = regexp(fstring,'STRUCTURE_ELEMENTS','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %f %f %f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);
strMsh.elements = out(:,2:4);
fprintf('>> Number of elements in the mesh: %d \n',length(strMsh.elements));

%% 8. Load the nodes on which homogeneous Dirichlet boundary conditions are applied
block = regexp(fstring,'STRUCTURE_DIRICHLET_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f','delimiter',' ','MultipleDelimsAsOne', 1);
    out{k} = horzcat(out{k}{:});
end
out = cell2mat(out);

% Filter out the actual DOFs
nDOFsPerNodeFromGiD = 3;
if strcmp(analysis.type,'PLANE_STRESS')||strcmp(analysis.type,'PLANE_STRAIN')
    nDOFsPerNode = 2;
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
    for j = 1:nDOFsPerNode
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

%% 9. Load the nodes on the Neumann boundary together with the load application information
block = regexp(fstring,'STRUCTURE_FORCE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %s %s','delimiter',' ','MultipleDelimsAsOne', 1);
end
out = out{1};
NBC.nodes = cell2mat(out(:,1));
outLoadType = out(:,2);
NBC.loadType = cell2mat(outLoadType{1});
outFctHandle = out(:,3);
NBC.fctHandle = cell2mat(outFctHandle{1});

%% 10. Get edge connectivity arrays for the Neumann edges
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
            [indexI,~] = find(nodeI == strMsh.elements);
            [indexJ,~] = find(nodeJ == strMsh.elements);
            
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

%% 11. Load the nodes on the IGA/FEM interface
block = regexp(fstring,'STRUCTURE_INTERFACE_NODES','split'); 
block(1) = [];
out = cell(size(block));
for k = 1:numel(block)
    out{k} = textscan(block{k},'%f %s %s','delimiter',' ','MultipleDelimsAsOne', 1);
end
out = out{1};
IBC.nodes = cell2mat(out(:,1));
% outLoadType = out(:,2);
% IBC.loadType = cell2mat(outLoadType{1});
% outFctHandle = out(:,3);
% IBC.fctHandle = cell2mat(outFctHandle{1});

%% 12. Get edge connectivity arrays for the IGA/FEM interface
fprintf('>> Interface boundary edges: %d \n',length(IBC.nodes)-1);

% Initialize the Interface boundary lines
IBC.lines = zeros(length(IBC.nodes)-1,3);

% Initialize line counter
counterLines = 1;

% Loop over each node pair
for i=1:length(IBC.nodes)
    for j=i:length(IBC.nodes)
        % If we are not in the same node
        if i~=j
            % Get the node index in the element array
            nodeI = IBC.nodes(i);
            nodeJ = IBC.nodes(j);
            
            
            % Find the element indices to which the nodes belong
            [indexI,~] = find(nodeI == strMsh.elements);
            [indexJ,~] = find(nodeJ == strMsh.elements);
            
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
                    IBC.lines(counterLines,:) = ...
                        [IBC.nodes(i) IBC.nodes(j) elementIndex];
                    
                    % Update counter
                    counterLines = counterLines + 1;
                end
            end
        end
    end
end

%% 13. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\nParsing took %.2d seconds \n\n',computationalTime);
    fprintf('_________________________Parsing Ended__________________________\n');
    fprintf('################################################################\n\n\n');
end

end