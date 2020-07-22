function flag = insidetri( A, B, C, p )
%INSIDETRI Checks if a point falls within the boundary of a triangle
%
% Inputs:
%       A   :
%       B   :
%       C   :
%       p   :
%
% Outputs:
%       flag: 
%
% Author: David Mostaza-Prieto
% The University of Manchester
% November 2012
%
%------------- BEGIN CODE --------------

v0 = C - A;
v1 = B - A;
v2 = p - A;

% Compute dot products
dot00 = dot(v0, v0);
dot01 = dot(v0, v1);
dot02 = dot(v0, v2);
dot11 = dot(v1, v1);
dot12 = dot(v1, v2);

% Compute barycentric coordinates
invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
u = (dot11 .* dot02 - dot01 .* dot12) .* invDenom;
v = (dot00 .* dot12 - dot01 .* dot02) .* invDenom;

% Check if point is inside triangle
flag = find((u>=0) & (v>=0) & (u+v<1));

%------------- END OF CODE --------------