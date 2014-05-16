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
function [u,noNewtonIter,hasConverged,FComplete,minElSize] = solve_FEMLinearSystem(uSaved, ...
    uDotSaved,uDDotSaved,mesh,analysis,F,bodyForces,parameters,u,uDot, ...
    computeLinearMatrices,DOFNumbering,freeDOFs,homDOFs,inhomDOFs, ...
    valuesInhomDOFs,strDynamics,nLinearAnalysis,int,outMsg)
%% Function documentation
%
% Returns the solution to a linear system which correspond to the
% classical finite element discretization of the underlying field.
%
%                    Input :
%                   uSaved : The discrete solution field of the previous 
%                            time step
%                uDotSaved : The time derivative of the discrete solution 
%                            field of the previous time step
%               uDDotSaved : The second order time derivative of the 
%                            discrete solution field of the previous time
%                            (dummy variable for this function)
%                     mesh : Nodes and elements of the underlying mesh
%                 analysis : On the analysis type
%                        F : The boundary force vector
%               bodyForces : The body force vector
%               parameters : The parameters of the physical field
%                        u : Initial guess for the primary field (just an 
%                            emtpy array must be given to this function)
%                     uDot : Initial guess for the time derivative of the
%                            primary field (dummy input for this function)
%    computeLinearMatrices : Function handle for the computation of the
%                            matrices and vectors nessecary for computing
%                            the solution vector
%             DOFNumbering : The global numbering of the DOFs in a
%                            3-dimensional array
%                 freeDOFs : The global numbering of the uncostrained DOFs
%                  homDOFs : The global numbering of the DOFs where
%                            homogeneous Dirichlet boundary conditions are
%                            applied
%                inhomDOFs : The global numbering of the DOFs where
%                            inhomogeneous Dirichlet boundary conditions 
%                            are applied
%          valuesInhomDOFs : The vector containing the magnitude of the
%                            prescribed values on the DOFs with
%                            inhomogeneous Dirichlet boundary conditions
%              strDynamics : Transient analysis parameters:
%                               .method : Time integration method
%                            .alphaBeta : Bossak parameter
%                                .gamma : Bossak parameter
%                                   .T0 : Start time of the simulation
%                                 .TEnd : End time of the simulation
%                                   .nT : Number of time steps
%                                   .dt : Time step (numeric or adaptive)
%          nLinearAnalysis : Nonlinear analysis parameters (dummy input 
%                            variable for this function)
%                               .scheme : The nonlinear solution scheme
%                            .tolerance : The residual tolerance
%                              .maxIter : The maximum number of nonlinear
%                                         iterations
%                      int : On the spatial integration
%                                      .type : 'default' or 'manual'
%                             .parentElement : 'tri', 'quad', usw.
%                                    .noXiGP : if .parentElement == 'quad'
%                                   .noEtaGP : if .parentElement == 'quad'
%                                      .noGP : if .parentElement == 'tri'
%                   outMsg : On printing information during analysis in the 
%                            command window
%
%                   Output :
%                        u : The converged discrete solution vector
%             noNewtonIter : No. of nonlinear iterations for convergence
%                            (dummy output for this function)
%             hasConverged : Flag on whether the nonlinear iterations have
%                            converged (for linear computations it is 
%                            returned always true)
%                FComplete : The complete force vector of the system
%                minElSize : The minimum element size over the isogeometric 
%                            mesh
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the linear matrices of the system
%
% 2. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
%
% 3. Solve the linear equation system
%
% 4. Re-assemble to the complete vector of unknowns
%
%% Function main body

%% 0. Read input

% Assign back dummy output variables
noNewtonIter = 'undefined';
hasConverged = true;

%% 1. Compute the linear matrices of the system
[K,RHS,minElSize] = computeLinearMatrices(u,uSaved,uDot,uDotSaved, ...
    DOFNumbering,mesh,analysis,F,bodyForces,strDynamics,parameters,int,outMsg);

%% 2. Update the right-hand side vector if inhomogeneous Dirichlet boundary conditions are encountered
if norm(valuesInhomDOFs) ~= 0
    RHS = RHS - K(:,inhomDOFs)*valuesInhomDOFs';
end

%% 3. Solve the linear equation system
uRed = sparse(K(freeDOFs,freeDOFs))\RHS(freeDOFs);

%% 4. Re-assemble to the complete vector of unknowns
u(freeDOFs) = uRed;
u(homDOFs) = 0;
u(inhomDOFs) = valuesInhomDOFs;

%% 5. Compute the complete force vector
FComplete = K*u;

end