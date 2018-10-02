function [idx, insert] = findRow(M, row)

% FINDROW Finds a given row in an ordered matrix.
%
%    idx = findRow(M, r) returns an index idx such that M(idx, :) == r, if it
%      is possible, otherwise idx == 0. The rows of the matrix M are supposed
%      to be lexicographically ordered. Running time is O(log n) where n is the
%      number of rows in the matrix.
%
%    [idx, insert] = findRow(M, r) moreover returns the index of the row where
%      r should be insered to keep M lexicographically sorted, if the row is not
%      in M. Otherwise, insert == 0.
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

if (isempty(M)),
    idx = 0;
    insert = 1;
else
    m = size(M, 2);
    low = 1;
    high = size(M, 1) + 1;
    
    % Bisection search
    while (low ~= high),
        mid = floor((low + high) / 2);
        
        for i=1:m+1,
            if (i == m + 1), % Line found
                idx = mid;
                insert = 0;
                return;
            end
            if row(i) < M(mid, i),
                high = mid;
                break;
            elseif row(i) > M(mid, i),
                low = mid + 1;
                break;
            end
        end
    end
    
    % Line not found
    idx = 0;
    insert = low;
end

