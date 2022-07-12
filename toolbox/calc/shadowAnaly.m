function [ shadPan ] = shadowAnaly( x, y, z, barC, delta, L_gw)
%SHADOWANALY Checks for mesh elements shadowed from a direction by others
%
% Inputs:
%       x       : X coordinates of the vertices of the N faces (3xN)
%       y       : Y coordinates of the vertices of the N faces (3xN)
%       z       : Z coordinates of the vertices of the N faces (3xN)
%       barC    : Barycentre of each mesh element (3xN)
%       delta   : Angle between flow and surface normal for each mesh element (1xN)
%       L_gw    : Rotation matrix from wind to geometric coordinate frame
%
% Outputs:
%       shadPan : Vector containing indices of the panels which are shadowed.
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
%
%--- Copyright notice ---%
% Copyright (C) 2021 The University of Manchester
% Written by David Mostaza Prieto,  Nicholas H. Crisp, Luciana Sinpetru and Sabrina Livadiotti
%
% This file is part of the ADBSat toolkit.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------

pAw = L_gw'*[x(1,:);y(1,:);z(1,:)];
pBw = L_gw'*[x(2,:);y(2,:);z(2,:)];
pCw = L_gw'*[x(3,:);y(3,:);z(3,:)];

barCw = L_gw'*barC;

xW = [pAw(1,:);pBw(1,:);pCw(1,:)];

xWmax = max(xW);
xWmin = min(xW);
 
indB = find(delta*180/pi>90.0001); % Backward-facing index (use 90.0001 to get a better split because panels at 90 degrees can't shadow)
indF = find(delta*180/pi<=90.0001); % Forward-facing index

minXwF = min(xWmin(indF));
maxXwB = max(xWmin(indB));

% Potential shadowed panels

indFPot = find(xWmin(indF)-maxXwB <10e-5); 
indBPot = find(xWmax(indB)-minXwF>10e-5);

% xWC = [pAw(1,indB(indBPot));pBw(1,indB(indBPot));pCw(1,indB(indBPot))];
yWC = [pAw(2,indB(indBPot));pBw(2,indB(indBPot));pCw(2,indB(indBPot))];
zWC = [pAw(3,indB(indBPot));pBw(3,indB(indBPot));pCw(3,indB(indBPot))];

shadPan = zeros(1,length(barCw));

tolB = 1e-5;
%tic
for i=1:length(indFPot)
    
    transYcord = yWC-barCw(2,indF(indFPot(i)));
%     tolY       = find(abs(transYcord<10e-6));
%     transYcord(tolY) = 0;
    
    yCoord   = abs(sum(sign(transYcord)));
    yC_chang = find(yCoord<3);
    
    if isempty(yC_chang)
        continue
    else
            transZcord = zWC(:,yC_chang)-barCw(3,indF(indFPot(i)));
%             tolZ       = find(abs(transZcord<10e-6));
%             transZcord(tolZ) = 0;
            
            zCoord   = abs(sum(sign(transZcord)));
            zC_chang = find(zCoord<3);
            
            if isempty(zC_chang)
                continue
            else

                % Check if [0,0] is inside the triangles
                
                p1 = [transYcord(1,yC_chang(zC_chang));transZcord(1,zC_chang)];
                p2 = [transYcord(2,yC_chang(zC_chang));transZcord(2,zC_chang)];
                p3 = [transYcord(3,yC_chang(zC_chang));transZcord(3,zC_chang)];
                
                Ftest = insidetri(p1,p2,p3,zeros(2,length(zC_chang)));
                
                if isempty(Ftest)

                    continue

                else
                 
                 if any((barCw(1,indB(indBPot((yC_chang(zC_chang(Ftest))))))-barCw(1,indF(indFPot(i))))> tolB)
                    shadPan(indF(indFPot(i))) = 1;
                 end
                 
                end
            end
    end
 
end

%toc
uj = find(shadPan>0);
shadPan = uj;

%------------- END OF CODE --------------
