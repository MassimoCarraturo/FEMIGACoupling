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
function [u,noNewtonIter,hasConverged,minElASize] = solve_IGANonlinearSystem(uSaved, ...
    uDotSaved,uDDotSaved,BSplinePatch,F,bodyForces,parameters,u,uDot, ...
    computeNonlinearMatrices,DOFNumbering,freeDOFs,homDOFs,inhomDOFs, ...
    valuesInhomDOFs,transientAnalysis,nLinearAnalysis,int,outMsg)
%% Function documentation
%
% Returns the solution to a nonlinear system which correspond to the
% isogeometric discretization of the underlying field.
%
%                    Input :
%                   uSaved : The discrete solution field of the previous 
%                            time step
%                uDotSaved : The time derivative of the discrete solution 
%                            field of the previous time step
%               uDDotSaved : The second order time derivative of the 
%                            discrete solution field of the previous time
%                            (dummy variable for this function)
%             BSplinePatch : The polynomial degrees, the knot vectors and 
%                            the Control Point coordinates/weights of the 
%                            B-Spline patch
%                        F : The boundary force vector
%               bodyForces : The body force vector
%               parameters : The parameters of the physical field
%                        u : Initial guess for the primary field
%                     uDot : Initial guess for the time derivative of the
%                            primary field
% computeNonlinearMatrices : Function handle for the computation of the
%                            matrices and vectors nessecary for the update
%                            of the discrete solution field
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
%        transientAnalysis : Transient analysis parameters:
%                               .method : Time integration method
%                            .alphaBeta : Bossak parameter
%                                .gamma : Bossak parameter
%                               .TStart : Start time of the simulation
%                                 .TEnd : End time of the simulation
%                                   .nT : Number of time steps
%                                   .dt : Time step (numeric or adaptive)
%          nLinearAnalysis : Nonlinear analysis parameters
%                               .method : The nonlinear solution scheme
%                                  .eps : The residual tolerance
%                              .maxIter : The maximum number of nonlinear
%                                         iterations
%                      int : On the spatial integration
%                                       .type : 'default' or 'manual'
%                                      .xiNGP : No. of GPs along xi-
%                                               direction for stiffness 
%                                               entries
%                                     .etaNGP : No. of GPs along eta-
%                                               direction for stiffness 
%                                               entries
%                               .xiNGPForLoad : No. of GPs along xi-
%                                               direction for load entries
%                              .etaNGPForLoad : No. of GPs along eta-
%                                               direction for load entries
%                                 .nGPForLoad : No. of GPs along boundary
%                   outMsg : On printing information during analysis in the 
%                            command window
%
%                   Output :
%                        u : The converged discrete solution vector
%             noNewtonIter : No. of nonlinear iterations for convergence
%             hasConverged : Flag on whether the nonlinear iterations hasve
%                            converged
%               minElASize : The minimum element area over the
%                            isogeometric mesh
%
% Function layout :
%
% 1. Loop over all the Newton-Rapson iterations
% ->
%    1i. Compute the stiffness matrix and the residual vector to the system
%
%   1ii. Compute the right-hand side (RHS) residual vector in equation
%
%  1iii. Check condition for convergence on the residual vector
%
%   1iv. Solve the linearized equation system
%
%    1v. Re-assemble to the complete vector of unknowns
% <-
%
% 2. Return the number of the nonlinear iterations needed for convergence
%
%% Function main body

%% 1. Loop over all the Newton-Rapson iterations
if strcmp(outMsg,'outputEnabled')
    msgPNR = sprintf('\t \t Looping over all the Newton iterations \n \t \t -------------------------------------- \n \n');
    fprintf(msgPNR);
end
for counterNewton = 1:nLinearAnalysis.maxIter
    %% 1i. Compute the stiffness matrix and the residual vector to the system
    [K,resVct,minElASize] = computeNonlinearMatrices(u,uSaved, ...
        uDot,uDotSaved,DOFNumbering,BSplinePatch,F,transientAnalysis, ...
        parameters,bodyForces,int);
    
    %% 1ii. Compute the right-hand side (RHS) residual vector in equation
    RHS = - resVct;
    if norm(valuesInhomDOFs(inhomDOFs)) ~= 0 && counterNewton == 1
        RHS = RHS + K(:,inhomDOFs)*valuesInhomDOFs(inhomDOFs);
    end
    
    %% 1iii. Check condition for convergence on the residual vector
    
    % Compute the norm of the residual vector over the free DOFs
    residualNorm = norm(resVct(freeDOFs));
    
    % Issue a message on the evolution of the residual vector
    if strcmp(outMsg,'outputEnabled')
        msgNR = sprintf('\t \t ||resVct|| = %d at nonlinear iteration No. = %d \n',residualNorm,counterNewton);
        fprintf(msgNR);
    end
    
    % Check the convergence of the Newton iterations
    if residualNorm<=nLinearAnalysis.eps && counterNewton ~= 1
        if strcmp(outMsg,'outputEnabled')
            msgANR = sprintf('\n \t \t Nonlinear iterations converged! \n \n \n');
            fprintf(msgANR);
        end
        hasConverged = true;
        break;
    end
    
    % If the Newton iterations do not converge after specified limit:
    if counterNewton == nLinearAnalysis.maxIter
        if strcmp(outMsg,'outputEnabled')
            if transientAnalysis.isAdaptive
                msgANR = sprintf('\n \t \t Nonlinear iterations did not converge, time step is adapted \n \n \n');
                fprintf(msgANR);
            else
                error('\n \t \t Nonlinear iterations did not converge for the fixed time step of %d',transientAnalysis.dt);
            end
        end

        % Flag on the convergence of the Newton iterations
        hasConverged = false;
    end
    
    %% 1iv. Solve the linearized equation system
    deltauRed = sparse(K(freeDOFs,freeDOFs))\RHS(freeDOFs);
    
    %% 1v. Re-assemble to the complete vector of unknowns
    u(freeDOFs) = u(freeDOFs) + deltauRed;
    if counterNewton == 1
        u(homDOFs) = 0;
        u(inhomDOFs) = valuesInhomDOFs(inhomDOFs);
    end
end

%% 2. Return the number of the nonlinear iterations needed for convergence
noNewtonIter = counterNewton;

end

