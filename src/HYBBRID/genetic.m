function [best_fit, best_arg, log] = genetic(fit_fun, breed_fun, data, params)

% GENETIC Genetic algorithm component of HYBBRID.
%
%    GENETIC(fit_fun, breed_fun, data, params) applies a genetic algorithm to
%      approximate the optimal OBB enclosing the set of points described in data
%      (one point per row). The term "optimal" corresponds to the minimal
%      fitness OBB, fit_fun being a handler to this fitness function. The
%      argument breed_fun is a handler to the function modifying the population
%      at each generation (crossover and mutation steps). More precisely,
%        - fit_fun(member, data, params) will be called to obtain the fitness
%          of a population member ;
%        - breed_fun(population, fitness, fit_fun, data, params) will be called
%          to obtain the next generation of the population, given the current
%          fitness values. In this case, both population and fitness are sorted
%          by increasing values of fitness.
%
%      The following parameters are used in this function:
%      params.g_popsize - population size
%      params.g_maxiter - maximum number of generations
%      params.g_tolval, params.g_toliter - stopping criterion: the algorithm
%                         is stopped if the relative improvement is < tolval
%                         during toliter iterations
%      params.g_verbose - display information at each generation if == 1
%
% output #1: The minimum fitness of a population member.
% output #2: The population member having the minimum fitness.
% output #3: A N-by-4 array containing for all iterations 1..N, the elapsed
%             time, the current minimum fitness value, the corresponding
%             population member, and the whole population.
%
% REMARKS:
%    A local optimization (localOptiRC) will be done when the current minimum
%    fitness value is improved. The initial population is generated at random.
%
% SEE ALSO:
%    HYBBRID, localOptiRC
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

% Randomly generate initial population
[n, p] = size(data);
m = params.g_popsize;
pop = cell(m, 1);

for i=1:m,
    pop{i} = cell(1, p+1);
    for j = 1:p+1,
        [q, r] = qr(rand(p, p)); %#ok<NASGU>
        pop{i}{j} = q;
    end
end

best_arg = pop{1};
best_fit = fit_fun(best_arg, data, params);
check_fit = Inf;
stop = 0;

% Evolution
starttime = clock;
log = cell(params.g_maxiter, 4);
for iter = 1:params.g_maxiter,
    % Compute fitness
    fitness = zeros(m, 2);
    for i = 1:m,
        [fitness(i, 1), details] = fit_fun(pop{i}, data, params);
        [foo, fitness(i, 2)] = min(details);
    end
    [fitsort, fitidx] = sortrows(fitness);
    old_fit = best_fit;
    % Check best fitness
    if (fitsort(1, 1) < check_fit),
        check_fit = fitsort(1, 1);
        [candidate, better_fit] = localOptiRC(data, pop{fitidx(1)}{fitsort(1, 2)});
        if (better_fit < best_fit),
            best_arg = candidate;
            best_fit = better_fit;
        end
    end

    % Create new population
    pop = breed_fun(pop(fitidx), fitsort, fit_fun, data, params);
    elapsedtime = etime(clock, starttime);
    if (params.g_verbose == 1),
        disp(sprintf('[%7.5g] Iter #%3d : objective value = %17.15g --- best value = %17.15g', elapsedtime, iter, fitsort(1), best_fit));
    end
    log{iter, 1} = elapsedtime;
    log{iter, 2} = best_fit;
    log{iter, 3} = best_arg;
    log{iter, 4} = pop;
    
    % Stopping criterion
    if (abs(best_fit - old_fit) < params.g_tolval*best_fit),
        stop = stop + 1;
        if (stop >= params.g_toliter),
            break;
        end
    else
        stop = 0;
    end
end

log = log(1:iter, :);

