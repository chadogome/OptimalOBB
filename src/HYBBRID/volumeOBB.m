function [volume, details] = volumeOBB(set, data, params)

% VOLUMEOBB Computation of the volumes of a set of OBBs enclosing a dataset.
%
%    VOLUMEOBB(set, data) computes the volumes of the OBBs enclosing the set of
%      points described in data (one point per row). The parameter set is a
%      cell array contatining rotation matrices associated to each OBB.
%
%    VOLUMEOBB(R, data) is the same as volumeOBB({R}, data), if R is a matrix.
%
% output #1: Volume of the smallest OBB among the given set.
% output #2: Volume of all OBBs in the set.
%
% REMARKS:
%    This function uses the MEX-function volumeAABBmex to compute the volume
%    of the OBB. If MEX is unavailable, it is possible to use the function
%    volumeAABB instead.
%
% SEE ALSO:
%    HYBBRID, volumeAABBmex
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

% Shortcut for sets with one rotation matrix
if (~iscell(set)),
    set = {set};
end

p = length(set);
details = zeros(1, p);
for i = 1:p,
%       details(i) = volumeAABB(data, set{i});        
        details(i) = volumeAABBmex(data, set{i}); % MEX version
end

volume = min(details);
