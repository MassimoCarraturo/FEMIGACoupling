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
function dR = computeNURBSBasisFunctionsAndDerivativesForCurve(knotSpanIndex,p,xi,Xi,CP,nDeriv)
%% Function documentation
%
% Returns an array containing the NURBS basis functions and their
% derivatives up to order nDeriv at the chosen parametric location. Source 
% reference:
%
% Les Piegl and Wayne Tiller, The NURBS Book. Springer-Verlag: Berlin 1995
% p. 72.
%
%         Input :
% knotSpanIndex : Knot span index
%             p : The polynomial degree of the NURBS curve
%            xi : The curve parameter where the basis functions are to be
%                 evaluated
%            Xi : The knot vector of the NURBS curve
%            CP : The set of Control Point coordinates and weights of the NURBS
%                 curve
%        nDeriv : The number of derivatives to be computed
%
%        Output :
%            dR : Array containing the NURBS basis functions and their
%                 derivatives up to order nDeriv. Element dR(i,j), i=1,...,
%                 p+1 and j=1,...,nDeriv+1 returns the (j+1)-th derivative 
%                 of the i-th non identically zero basis function at the
%                 knot span with index knotSpanIndex
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the BSpline basis functions and their derivatives at xi
%
% 2. Compute the denominator function and its derivatives at xi
%
% 3. Loop over all the basis functions
%
%    3i. Loop over all the derivatives
%
%        3i.2. Compute the Control point index
%
%        3i.3. Compute the product of the derivatives of the basis functions with the Control Point weights
%
%        3i.4. Loop over all the involved derivatives
%
%        3i.5. Divide by the denominator function
%
%% Function main body

%% 0. Read input

% Initialize the output array
dR = zeros(p+1,nDeriv+1);

%% 1. Compute the BSpline basis functions and their derivatives at xi
dN = computeBSplineBasisFunctionsAndDerivativesForCurve(knotSpanIndex,p,xi,Xi,nDeriv);

%% 2. Compute the denominator function and its derivatives at xi
dF = computeDenominatorFunctionAndDerivativesForCurve(dN,knotSpanIndex,p,CP,nDeriv);

%% 3. Loop over all the basis functions
for i=0:p
    %% 3i. Loop over all the derivatives
    for j=0:nDeriv
        %% 3i.2. Compute the Control point index
        index = knotSpanIndex-p+i;
        
        %% 3i.3. Compute the product of the derivatives of the basis functions with the Control Point weights
        v = dN(i+1,j+1)*CP(index,4);
        
        %% 3i.4. Loop over all the involved derivatives
        for k=1:j
            v = v - nchoosek(j,k)*dF(k+1,1)*dR(i+1,j-k+1);
        end
        
        %% 3i.5. Divide by the denominator function
        dR(i+1,j+1) = v/dF(1,1);
    end
end

end