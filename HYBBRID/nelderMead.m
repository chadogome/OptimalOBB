function simplex = nelderMead(simplex, simplex_fun, data, params)

% NELDERMEAD Nelder-Mead simplex algorithm on SO(3, R)
%
%    NELDERMEAD(simplex, simplex_fun, data, params) applies a Nelder-Mead
%      algorithm with simplex as starting point. A handler to the function to
%      minimize is given by simplex_fun. More precisely, simplex_fun(simplex,
%      data, params) will be called to obtain the value of a simplex. The only
%      parameter used by this function is params.nm_maxiter, which is the
%      number of iterations that will be done. Default value is 20.
%
% output #1: The simplex after params.nm_maxiter iterations of Nelder-Mead.
%
% SEE ALSO:
%    HYBBRID
%
% AUTHORS:
%    Chia-Tche Chang <chia-tche.chang@uclouvain.be>
%    Bastien Gorissen <bastien.gorissen@cenaero.be>
%    Samuel Melchior <samuel.melchior@uclouvain.be>
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

p = length(simplex);
N = p-1;
if (nargin < 4), params = struct; end
if (~isfield(params, 'nm_maxiter')),   params.nm_maxiter = 20; end

% Standard values for Nelder-Mead simplex algorithm
rho = 1/2;
sigma = 1/2;
Rg = simplex{1};

for iter = 1:params.nm_maxiter,
    % Step 1 : Reordering
    [foo, simplex_val] = simplex_fun(simplex, data, params);
    [fi, idx] = sort(simplex_val);
    simplex = simplex(idx);

    % Step 2 : Computation of the center of gravity
    Rg = karcher(Rg, simplex(1:N));
    Rr = Rg*simplex{p}'*Rg;
    fr = simplex_fun(Rr, data, params);
    if (fr < fi(N)),
        if (fr >= fi(1)),
            % Step 3 : Reflection
            simplex{p} = Rr;
        else
            % Step 4 : Expansion
            Re = Rg*simplex{p}'*Rr;
            fe = simplex_fun(Re, data, params);
            if (fe < fr),
                simplex{p} = Re;
            else
                simplex{p} = Rr;
            end
        end
    else
        % Step 5 : Contraction
        Rc = affine(Rg, simplex{p}, rho, Rg, -rho, simplex{p});
        fc = simplex_fun(Rc, data, params);
        if (fc <= fi(p)),
            simplex{p} = Rc;
        else
            % Step 6 : Reduction
            for i = 2:p,
                simplex{i} = affine(Rg, simplex{1}, sigma, simplex{i}, -sigma, simplex{1});
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kmean = karcher(P, X)

% Karcher mean of Xi's with respect to the point P

kmean = 0;
for i = 1:length(X),
    kmean = kmean + X{i};
end
kmean = kmean/length(X);
[q, r] = qr(kmean);
kmean = q/det(q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = affine(foot, base, factd, dest, facto, orig)

% Affine combination

res = zeros(size(foot));
if (~isempty(base)), res = res + base; end
if (~isempty(dest)), res = res + factd*dest; end
if (~isempty(orig)), res = res + facto*orig; end
[q, r] = qr(res);
res = q/det(q);

