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
function dR = computeIGABasisFunctionsAndDerivativesForSurface(xiKnotSpan,p,xi,Xi,etaKnotSpan,q,eta,Eta,CP,isNURBS,mixedDerivOrder)
%% Function documentation
%
% Returns the IGA basis functions (B-Spline or NURBS) and their derivatives 
% for the parametrization of a surface.
%
% The basis functions array dR(Rxi,Reta,derivs) is sorted as:
%
% dR(:,:) = [dR(1,1) dR(2,1) ... dR(xiNoBasisFct,1) dR(1,2) dR(2,2) ...
%            dR(xiNoBasisFct,2) ... dR(1,etaNoBasisFct) dR(2,etaNoBasisFct) 
%            ... dR(xiNoBasisFct,etaNoBasisFct)]'
%
% where dR denotes the derivatives. The derivatives themselves are sorted
% into the second dimension of the array which given the maximum order of 
% the derivatives n = mixedDerivOrder is given by:
%
% dR(derivs) = [R dR/dxi d^2R/dxi^2 ... d^nR/dxi^n 
%               dR/deta d^2R/deta*dxi ... d^(n-1)R/deta*dxi^(n-1) 
%               dR^2/deta^2 d^3R/deta^2*dxi ... d^(n-2)R/deta^2*dxi^(n-2)
%               ...                         ... ...
%               dR^(n-1)/deta^(n-1) d^(n-1)R/deta^(n-1)*dxi 
%               dR^n/deta^n]
%
% To get the i-th derivative with respect to xi-direction and j-th
% derivative with respect to eta-direction use the index 
%
% derivIndex = (mixedDerivOrder - j)*(mixedDerivOrder - j + 1) + i
%
% and to get the k-th basis function with its derivatives themselves use
% the index
%
% derivIndex = 
%            computeIndexForBSplineBasisFunctionsAndDerivatives(nDeriv,i,j)
%
% namely in total dR(basisFncIndex,derivIndex)
%
%           Input :
%      xiKnotSpan : The knot span index in xi-direction
%     etaKnotSpan : The knot span index in eta-direction
%               p : The polynomial degree in xi-direction
%               q : The polynomial degree in eta-direction
%          Xi,Eta : The knot vectors in xi-,eta- direction
%              xi : The surface parameter in xi-direction
%             eta : The surface parameter in eta-direction
%              CP : The set of the control point coordinates and weights
%         isNURBS : Flag on wheter the basis is a NURBS or a B-Spline
% mixedDerivOrder : The maximum order for the mixed derivative
%
%          Output :
%              dR : The array containing the B-Spline/NURBS basis functions 
%                   and their derivatives sorted as explained above
%
%% Function main body

if isNURBS
    dR = computeNURBSBasisFunctionsAndDerivativesForSurface(xiKnotSpan,p,xi,Xi,etaKnotSpan,q,eta,Eta,CP,mixedDerivOrder);
else
    dR = computeBSplineBasisFunctionsAndDerivativesForSurface(xiKnotSpan,p,xi,Xi,etaKnotSpan,q,eta,Eta,mixedDerivOrder);
end

end