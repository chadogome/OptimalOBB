function [Ropt, Vopt] = localOptiRC(data, R)

% LOCALOPTIRC Local optimization of an OBB in the three directions.
%
%    LOCALOPTIRC(data, R) computes the OBB of minimum volume in the neigh-
%      borhood of the bounding box oriented w.r.t. the rotation matrix R. The
%      neighborhood consists of all OBBs that may be obtained by applying the
%      2D Rotating Calipers technique on one of the three axes defined by R.
%
% output #1: Rotation matrix associated to the optimal OBB in the neighborhood.
% output #2: Volume of the optimal OBB in the neighborhood.
%
% SEE ALSO:
%    HYBBRID, rotatingCalipers
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

% Preprocessing step
idx = convhulln(data);
data = data(unique(idx), :);
rotatedData = data*R';

% Apply the rotating calipers on the three axes
rangeData = max(rotatedData) - min(rotatedData);
[volx, anglex] = rotatingCalipers(rotatedData(:, [2, 3]));
[voly, angley] = rotatingCalipers(rotatedData(:, [1, 3]));
[volz, anglez] = rotatingCalipers(rotatedData(:, [1, 2]));
volx = volx*rangeData(1);
voly = voly*rangeData(2);
volz = volz*rangeData(3);

% Return the best OBB
[Vopt, idx] = min([volx, voly, volz]);
switch(idx),
    case 1,
        Ropt = [1 0 0 ; 0 cos(anglex) sin(anglex) ; 0 -sin(anglex) cos(anglex)]*R;
    case 2,
        Ropt = [cos(angley) 0 -sin(angley) ; 0 1 0 ; sin(angley) 0 cos(angley)]*R;
    case 3,
        Ropt = [cos(anglez) sin(anglez) 0 ; -sin(anglez) cos(anglez) 0 ; 0 0 1]*R;
end

