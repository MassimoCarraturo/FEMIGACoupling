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
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeOutputFEMPlateInMembraneActionToVTK(strMsh,nodalDisplacement,epsilon,sigma,caseName,pathToOutput,title)
%% Function documentation 
%
% Writes out the results of a classical FE plate in membrane action
% analysis into a VTK file to be read by paraview.
%
%             Input :
%            strMsh : Nodes and elements in the mesh
% nodalDisplacement : The displacement field arranged into a 2D array 
%                     [[x-comp ycomp],noNodes]
%           epsilon : The strain field [3,[epsilonXX epsilonYY epsilonXY]]
%             sigma : The stress field [3,[sigmaXX sigmaYY sigmaXY]]
%          caseName : The name of the case in the inputGiD case folder
%      pathToOutput : Path to the output file
%    outputFilename : The name of the output file
%
%            Output :
%                     Write results into file
%
%
%% Function main body

%% 0. Read input

%  Number of nodes in the mesh
[noNodes,~] = size (strMsh.nodes);

% Number of elements in the mesh
[noElements,elementOrder] = size(strMsh.elements);

%  Open the output file
if isempty(caseName)
    caseName = 'ns2d_fem.vtk';
end

% outputUnit = fopen(strcat(pathToOutput,outputFilename),'w' );
outputUnitColorPlots = fopen(strcat(pathToOutput,caseName,'/',caseName,'_contourPlots.vtk'),'w');
outputUnitDeformed = fopen(strcat(pathToOutput,caseName,'/',caseName,'_deformation.vtk'),'w');

% Transpose the nodal coordinates array
XYZ = strMsh.nodes(1:noNodes,:)';

if ( elementOrder == 6 )
	fprintf ( 1, '\n' );
    fprintf ( 1, 'TWO_TO_VTK:\n' );
    fprintf ( 1, '  The input data uses quadratic elements.\n' );
    fprintf ( 1, '  The output data will use linear elements.\n' );
end

% Re-arrange the element numbering to start from zero
elements = zeros(elementOrder,noElements);
elements(1:elementOrder,1:noElements) = strMsh.elements(1:noElements,1:elementOrder)' - 1;

%% 1. Write out the data for the color plots

% Write out the preamble
fprintf(outputUnitColorPlots, '# vtk DataFile Version 2.0\n');
fprintf(outputUnitColorPlots, '%s\n',title);
fprintf(outputUnitColorPlots, 'ASCII\n');
fprintf(outputUnitColorPlots, '\n');
fprintf(outputUnitColorPlots, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(outputUnitColorPlots, 'POINTS %d double\n', noNodes);

fprintf(outputUnitDeformed, '# vtk DataFile Version 2.0\n');
fprintf(outputUnitDeformed, '%s\n',title);
fprintf(outputUnitDeformed, 'ASCII\n');
fprintf(outputUnitDeformed, '\n');
fprintf(outputUnitDeformed, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(outputUnitDeformed, 'POINTS %d double\n', noNodes);

% Write out the nodal coordinates
for nodeID = 1:noNodes
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n', XYZ(:,nodeID));
    fprintf(outputUnitDeformed,'  %f  %f  %f\n', XYZ(:,nodeID) + nodalDisplacement(:,nodeID));
end

% Note that CELL_SIZE uses ELEMENT_ORDER+1 because the order of each 
% element is included as a data item.
cellSize = noElements*(elementOrder + 1);

% Output the element connectivities to the nodes
fprintf(outputUnitColorPlots,'\n' );
fprintf(outputUnitColorPlots,'CELLS  %d  %d\n',noElements,cellSize);

fprintf(outputUnitDeformed,'\n' );
fprintf(outputUnitDeformed,'CELLS  %d  %d\n',noElements,cellSize);

% Loop over all the elements in the mesh
for elementID = 1:noElements
    % Write out the element order
    fprintf (outputUnitColorPlots,'  %d',elementOrder);
    fprintf (outputUnitDeformed,'  %d',elementOrder);

    % Loop over all the polynomial orders
    for order = 1:elementOrder
        % Write out the polynomial order
        fprintf(outputUnitColorPlots,'  %d',elements(order,elementID));
        fprintf(outputUnitDeformed,'  %d',elements(order,elementID));
    end
    
    % Change line
    fprintf (outputUnitColorPlots,'\n' );
    fprintf (outputUnitDeformed,'\n' );
end

% VTK has a cell type 22 for quadratic triangles.  However, we
% are going to strip the data down to linear triangles for now,
% which is cell type 5.

fprintf (outputUnitColorPlots,'\n' );
fprintf (outputUnitColorPlots,'CELL_TYPES %d\n',noElements);

fprintf (outputUnitDeformed,'\n' );
fprintf (outputUnitDeformed,'CELL_TYPES %d\n',noElements);

% Loop over all the elements and write out the nodal coordinates and the
% element connectivities according to the element order
if elementOrder == 3
    for elementID = 1 : noElements
        fprintf (outputUnitColorPlots,'5\n');
        fprintf (outputUnitDeformed,'5\n');
    end
elseif elementOrder == 6
    for elementID = 1:noElements
        fprintf (outputUnitColorPlots,'22\n');
        fprintf (outputUnitDeformed,'22\n');
    end
end

% Write out the strain tensor
fprintf(outputUnitColorPlots,'\n' );
fprintf(outputUnitColorPlots,'CELL_DATA %d\n',noElements);
fprintf(outputUnitColorPlots,'TENSORS strain double\n');
for elementID = 1 : noElements
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n',epsilon(1,elementID),epsilon(3,elementID),0.0);
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n',epsilon(3,elementID),epsilon(2,elementID),0.0);
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n',0.0,0.0,0.0);
    if elementID~=noElements
        fprintf(outputUnitColorPlots,'\n');
    end
end

% Write out the displacement vector field at each node
fprintf(outputUnitColorPlots,'POINT_DATA %d\n',noNodes);
fprintf(outputUnitColorPlots,'VECTORS displacements double\n');
for nodeID = 1:noNodes
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n',nodalDisplacement(:,nodeID));
end

% Write out the stress tensor
fprintf(outputUnitColorPlots,'\n' );
fprintf(outputUnitColorPlots,'CELL_DATA %d\n',noElements);
fprintf(outputUnitColorPlots,'TENSORS stress double\n');
for elementID = 1 : noElements
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n',sigma(1,elementID),sigma(3,elementID),0.0);
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n',sigma(3,elementID),sigma(2,elementID),0.0);
    fprintf(outputUnitColorPlots,'  %f  %f  %f\n\n',0.0,0.0,0.0);
    if elementID~=noElements
        fprintf(outputUnitColorPlots,'\n');
    end
end

%% 2. Close files
fclose(outputUnitColorPlots);
fclose(outputUnitDeformed);

return

end
