function flag = insidetri( A, B, C, p )
% Checks if a point falls within the boundary of a triangle.
% Assumes the triangles have been projected onto a 2-dimensional plane,
% and therefore all points have [x,y] coordinates only.
%
% Inputs:
%       A   : 2xN matrix, each column is the [x,y] coordinates of a point of the triangle
%       B   : 2xN matrix, each column is the [x,y] coordinates of the second point of the triangle
%       C   : 2xN matrix, each column is the [x,y] coordinates of the third point of the triangle
%       p   : 2xN matrix, each column is the [x,y] coordinates of a point to be checked for inclusion in the triangle [A,B,C]
%
% Outputs:
%       flag: Vector containing the index of any points which fell inside their respective triangles.
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
