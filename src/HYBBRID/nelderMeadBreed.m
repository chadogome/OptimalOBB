function newpop = nelderMeadBreed(pop, fitness, fit_fun, data, params)

% NELDERMEADBREED Genetic population breed using Nelder-Mead (HYBBRID).
%
%    NELDERMEADBREED(pop, fitness, fit_fun, data, params) generates a new
%      population given the current population described by the cell array pop,
%      the fitness value of each population member (double array fitness) and a
%      handler to this fitness function fit_fun. More precisely,
%      fit_fun(member, data, params) will be called to obtain the fitness of a
%      population member. Each population member is expected to be a simplex of
%      rotation matrices, and the argument data represents the associated set of
%      points (one point per row). A simplex is represented by a cell array
%      containing several double arrays. Details about the population generation
%      can be found in the reference. The only parameter used in this function
%      is params.g_randmut, which corresponds to probability to have random
%      mutations occuring at any given generation (default value is 0). If
%      necessary, such random mutations are applied in the following way: let m
%      be the population size, then a population member is randomly chose and
%      replaced by a random simplex. This is repeated m times.
%
% output #1: The new population of simplices.
%
% SEE ALSO:
%    HYBBRID, nelderMead
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

if (nargin < 5), params = struct; end
if (~isfield(params, 'g_randmut')),   params.g_randmut = 0.0; end

m = length(pop);
p = length(pop{1});
mmf = floor(m/2);
mmc = ceil(m/2);

% Creation of the four groups
idx = ceil(mmc*rand(mmf, 1));
pop1 = pop(idx);
fit1 = fitness(idx);
idx = ceil(mmc*rand(mmf, 1));
pop2 = pop(idx);
fit2 = fitness(idx);
idx = ceil(mmf*rand(mmc, 1));
pop3 = pop(idx);
fit3 = fitness(idx);
idx = ceil(mmf*rand(mmc, 1));
pop4 = pop(idx);
fit4 = fitness(idx);

% Crossover I: pop1 x pop2
newpop1 = cell(mmf, 1);
cutoff = 0.5 + 0.1*(fit1 <= fit2) - 0.1*(fit1 >= fit2);
for i=1:mmf,
    newpop1{i} = cell(1, p);
    for j=1:p,
        if (rand < cutoff(i)),
            newpop1{i}{j} = pop1{i}{j};
        else
            newpop1{i}{j} = pop2{i}{j};
        end
    end
end

% Crossover II: pop3 x pop4
newpop2 = cell(mmc, 1);
cutoff = 0.5 + 0.1*(fit3 <= fit4) - 0.1*(fit3 >= fit4);
for i=1:mmc,
    newpop2{i} = cell(1, p);
    for j=1:p,
        newpop2{i}{j} = affine(pop3{i}{j}, [], cutoff(i), pop3{i}{j}, 1-cutoff(i), pop4{i}{j});
    end
end

% Nelder-Mead mutation
newpop = [newpop1 ; newpop2];
for i = 1:m,
    newpop{i} = nelderMead(newpop{i}, fit_fun, data, params);
end

% Random mutations
if rand < params.g_randmut,
    for ii = 1:m,
        i = ceil(m*rand);
        newpop{i} = cell(1, p);
        for j = 1:p,
            [q, r] = qr(rand(p-1, p-1)); %#ok<NASGU>
            newpop{i}{j} = q;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = affine(foot, base, factd, dest, facto, orig)

res = zeros(size(foot));

if (~isempty(base)), res = res + base; end
if (~isempty(dest)), res = res + factd*dest; end
if (~isempty(orig)), res = res + facto*orig; end

[q, r] = qr(res);
res = q/det(q);

