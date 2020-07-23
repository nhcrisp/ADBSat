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
%       shadPan : Flag for shadowed mesh elements (1xN)
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
%
%------------- BEGIN CODE --------------

pAw = L_gw'*[x(1,:);y(1,:);z(1,:)];
pBw = L_gw'*[x(2,:);y(2,:);z(2,:)];
pCw = L_gw'*[x(3,:);y(3,:);z(3,:)];

barCw = L_gw'*barC;

xW = [pAw(1,:);pBw(1,:);pCw(1,:)];

xWmax = max(xW);
xWmin = min(xW);

indB = find(delta*180/pi>(90.01)); % Backward-facing index
indF = find(delta*180/pi<=(90.01)); % Forward-facing index

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
                 
                 if any((barCw(1,indB(indBPot((yC_chang(zC_chang(Ftest))))))-barCw(1,indF(indFPot(i))))> tolB);
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
