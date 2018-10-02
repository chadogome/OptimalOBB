function vol = volumeAABB(data, R, initCache)

% VOLUMEAABB Computation of the volume of the minimal AABB enclosing a dataset.
%
%    VOLUMEAABB(data) computes the volume of the minimal axis-aligned bounding
%      box enclosing the set of points described in data (one point per row).
%
%    VOLUMEAABB(data, R) is the same as volumeAABB(data*R').
%
% output #1: Volume of the AABB enclosing data.
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

if (nargin == 2),
    data = data*R';
end
vol = prod(max(data) - min(data));

