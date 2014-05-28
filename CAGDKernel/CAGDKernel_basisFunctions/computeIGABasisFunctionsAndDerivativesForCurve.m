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
function dR = computeIGABasisFunctionsAndDerivativesForCurve(knotSpanIndex,p,xi,Xi,CP,nDeriv,isNURBS)
%% Function documentation
%
% Returns an array containing the B-Spline/NURBS basis functions and their
% derivatives up to order nDeriv at the chosen parametric location.
%
%         Input :
% knotSpanIndex : Knot span index
%             p : The polynomial degree of the B-Spline curve
%            xi : The curve parameter where the basis functions are to be
%                 evaluated
%            Xi : The knot vector of the B-Spline curve
%            CP : The set of Control Point coordinates and weights of the 
%                 B-Spline curve
%        nDeriv : The number of derivatives to be computed
%       isNURBS : Flag on the whether the basis is a B-Spline or a NURBS
%
%        Output :
%            dR : Array containing the B-Spline/NURBS basis functions and 
%                 their derivatives up to order nDeriv. Element dR(i,j), 
%                 i=1,...,p+1 and j=1,...,nDeriv+1 returns the (j+1)-th 
%                 derivative of the i-th non identically zero basis 
%                 function at the knot span with index knotSpanIndex
%
%% Function main body

if isNURBS
    dR = computeNURBSBasisFunctionsAndDerivativesForCurve(knotSpanIndex,p,xi,Xi,CP,nDeriv);
else
    dR = computeBSplineBasisFunctionsAndDerivativesForCurve(knotSpanIndex,p,xi,Xi,nDeriv);
end

end

