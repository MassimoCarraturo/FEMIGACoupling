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
function [xs,ys,zs] = createSupports2D(CP,rb)
%% Function documentation
%
% returns coordinates for the triangles at the support locations
%
%   Input : 
%      CP : Control Point coordinates
%      rd : Array containing information on the supports
%
%  Output :
%      xs : x-coordinates of the support triangle vertices
%      ys : y-coordinates of the support triangle vertices
%      zs : z-coordinates of the support triangle vertices
%
%% Function main body

% Total number of Control Points
nu = length(CP(:,1,1));

% scaling factors for the support triangles
up = max(max(max(max(CP))));
lo = min(min(min(min(CP))));
% Average the factor with respect to the maximum and minimum values 
fac = (up-lo)/5;

% Initialize the output arrays
xs = zeros(length(rb),4);
ys = zeros(length(rb),4);
zs = zeros(length(rb),4);

for k = 1:length(rb)
    % Get the corresponding Control Point number p and indices CP(i,j)
    h = rb(k)/2;
    p = ceil(h);
    j = ceil(p/nu);
    i = p-(j-1)*nu;
    %(rb is odd -> horizontal support)
    if (p~=h)   
        xs(k,1)=CP(i,j,1);
        xs(k,2)=CP(i,j,1)-0.1732*fac;
        xs(k,3)=CP(i,j,1)-0.1732*fac;
        xs(k,4)=xs(k,1);
        ys(k,1)=CP(i,j,2);
        ys(k,2)=CP(i,j,2)+0.1*fac;
        ys(k,3)=CP(i,j,2)-0.1*fac;
        ys(k,4)=ys(k,1);
        zs(k,1:4)=CP(i,j,3);
    %(rb is even -> vertical support)
    else        
        xs(k,1)=CP(i,j,1);
        xs(k,2)=CP(i,j,1)-0.1*fac;
        xs(k,3)=CP(i,j,1)+0.1*fac;
        xs(k,4)=xs(k,1);
        ys(k,1)=CP(i,j,2);
        ys(k,2)=CP(i,j,2)-0.1732*fac;
        ys(k,3)=CP(i,j,2)-0.1732*fac;
        ys(k,4)=ys(k,1);
        zs(k,1:4)=CP(i,j,3);
    end
end

end