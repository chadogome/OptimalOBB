function [volume, angle] = rotatingCalipers(data)

% ROTATINGCALIPERS Minimal bounding rectangle for 2D datasets.
%
%    ROTATINGCALIPERS(data) computes the minimal oriented bounding rectangle
%      enclosing the points in R^2 described in data (one point per row).
%
% output #1: Volume of the minimal oriented bounding rectangle.
% output #2: Orientation of the minimal oriented bounding rectangle, i.e.,
%             angle between one edge of the box and the x-axis, in radians.
%
% SEE ALSO:
%    HYBBRID
%
% AUTHORS:
%    Chia-Tche Chang <cchang.uclouvain@gmail.com>
%    Bastien Gorissen <bastien@panopticgame.com>
%    Samuel Melchior <samuel.melchior@epfl.ch>
%
% REFERENCES:
%    - C.-T.Chang, B.Gorissen, S.Melchior, 
%      "Fast oriented bounding box optimization on the rotation group SO(3, R)",
%      submitted to ACM Transactions on Graphics, 2011.
%    - G. Toussaint, "Solving geometric problems with the rotating calipers",
%      Proc. IEEE MELECON '83, 10-02, 1983.

% This file is part of HYBBRID.
% Copyright © 2010 Chia-Tche Chang, Bastien Gorissen, Samuel Melchior
%
% HYBBRID is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% HYBBRID is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with HYBBRID. If not, see <http://www.gnu.org/licenses/>.
% Computation of the convex hull
convHull = convhull(data(:, 1), data(:, 2));
convHull = convHull(1:end-1, :);
hullData = data(convHull, :);
nbVertices = size(hullData, 1);
[minxy, argminxy] = min(hullData);
[maxxy, argmaxxy] = max(hullData);

% Initialization
Theta = 0;
Right = argmaxxy(1);
Up = argmaxxy(2);
Left = argminxy(1);
Down = argminxy(2);
NextRight = mod(Right, nbVertices) + 1;
NextUp = mod(Up, nbVertices) + 1;
NextLeft = mod(Left, nbVertices) + 1;
NextDown = mod(Down, nbVertices) + 1;
DeltaRight = hullData(NextRight, :) - hullData(Right, :);
DeltaUp = hullData(NextUp, :) - hullData(Up, :);
DeltaLeft = hullData(NextLeft, :) - hullData(Left, :);
DeltaDown = hullData(NextDown, :) - hullData(Down, :);
AngleRight = acos(DeltaRight(2)/norm(DeltaRight));
AngleUp = acos(-DeltaUp(1)/norm(DeltaUp));
AngleLeft = acos(-DeltaLeft(2)/norm(DeltaLeft));
AngleDown = acos(DeltaDown(1)/norm(DeltaDown));

volume = Inf;
angle = NaN;

% Rotate the calipers
while (2*Theta < pi),
    rotatedData = hullData*[cos(Theta), -sin(Theta) ; sin(Theta), cos(Theta)];
    [minxy, argminxy] = min(rotatedData);
    [maxxy, argmaxxy] = max(rotatedData);
    curVolume = (maxxy(1) - minxy(1))*(maxxy(2) - minxy(2));
    if (curVolume < volume),
        volume = curVolume;
        angle = Theta;
    end
    
    [Theta, argminth] = min([AngleRight, AngleUp, AngleLeft, AngleDown]);
    switch(argminth),
        case 1,
            Right = NextRight;
            NextRight = mod(Right, nbVertices) + 1;
            DeltaRight = hullData(NextRight, :) - hullData(Right, :);
            AngleRight = acos(DeltaRight(2)/norm(DeltaRight));
        case 2,
            Up = NextUp;
            NextUp = mod(Up, nbVertices) + 1;
            DeltaUp = hullData(NextUp, :) - hullData(Up, :);
            AngleUp = acos(-DeltaUp(1)/norm(DeltaUp));
        case 3,
            Left = NextLeft;
            NextLeft = mod(Left, nbVertices) + 1;
            DeltaLeft = hullData(NextLeft, :) - hullData(Left, :);
            AngleLeft = acos(-DeltaLeft(2)/norm(DeltaLeft));
        case 4,
            Down = NextDown;
            NextDown = mod(Down, nbVertices) + 1;
            DeltaDown = hullData(NextDown, :) - hullData(Down, :);
            AngleDown = acos(DeltaDown(1)/norm(DeltaDown));
    end
end

