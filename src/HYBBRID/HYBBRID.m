function [best, bestvalue, log, data] = HYBBRID(n, params)

% HYBBRID HYbrid Bounding Box Rotation IDentification.
%
%    HYBBRID(n) approximates the minimum-volume Oriented Bounding Box (OBB)
%      enclosing a random set of n points in the 3D unit cube (default n = 50).
%
%    HYBBRID(data) approximates the minimum-volume OBB enclosing the set of
%      points described in data (one point per row in the matrix).
%
%    HYBBRID(..., params) does the same as HYBBRID(...) with customized
%      parameters (default values = {1, 1, 30, 100, 1e-2, 5, 1, 0.0, 20}):
%      params.test_repeat  - repeat the execution N times
%                            (benchmark mode if N != 1)
%      params.opt_convhull - preprocess data by extracting convex hull if == 1
%      params.g_popsize    - population size of the genetic algorithm
%      params.g_maxiter    - maximum number of generations of the genetic
%                            algorithm
%      params.g_tolval, params.g_toliter - stopping criterion: the genetic
%                            algorithm is stopped if the relative improvement is
%                            < tolval during toliter iterations
%      params.g_verbose    - display information at each generation if == 1
%      params.g_randmut    - apply random mutations with given probability at
%                            each generation
%      params.nm_maxiter   - number of Nelder-Mead iterations at each generation
%
% output #1: Rotation matrix associated to the orientation of the minimum OBB
%            found by the algorithm.
% output #2: Volume of the minimum OBB found by the algorithm.
% output #3: In benchmark mode, this is a double array containing, for each
%            five iterations, the min-max-mean values of the currently elapsed
%            time, the current best value and the number of runs having reached
%            this iteration, in a row.
%            Three additional rows contain the following information:
%            - The min-max-mean of the total elapsed time and the min-max-min-
%              std of the minimum volume found by the algorithm.
%            - The {5-10-25-50-75-90-95}th percentiles of the total elapsed
%              time.
%            - The {5-10-25-50-75-90-95}th percentiles of the minimum volume
%              found by the algorithm.
%            In normal mode, the current elapsed time, current best volume
%            and current best rotation matrix (for each iteration) are returned
%            in a cell array.
% output #4: Matrix representing the dataset (after preprocessing if enabled).
%
% SEE ALSO:
%    genetic, nelderMeadBreed
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

if (nargin == 0),
    n = 50;
elseif (isvector(n) == 1),
    n = sum(n);
end

if (isscalar(n)),
    data = rand(n, 3);
else
    data = n;
end

if (nargin < 2),
    params = struct;
end

% Default parameters
if (~isfield(params, 'test_repeat')),   params.test_repeat = 1; end
if (~isfield(params, 'opt_convhull')),  params.opt_convhull = 1; end
if (~isfield(params, 'g_popsize')),     params.g_popsize = 30; end
if (~isfield(params, 'g_maxiter')),     params.g_maxiter = 100; end
if (~isfield(params, 'g_tolval')),      params.g_tolval = 1e-2; end
if (~isfield(params, 'g_toliter')),     params.g_toliter = 5; end
if (~isfield(params, 'g_verbose')),     params.g_verbose = (params.test_repeat == 1); end
if (~isfield(params, 'g_randmut')),     params.g_randmut = 0.0; end
if (~isfield(params, 'nm_maxiter')),    params.nm_maxiter = 20; end

starttime = clock;
% Preprocessing with convex hull (default)
if (params.opt_convhull == 1),
    data = data(unique(convhulln(data)), :);
end

if (params.test_repeat == 1),
    % Genetic Nelder-Mead
    [bestvalue, best, log] = genetic(@volumeOBB, @nelderMeadBreed, data, params);
    disp(sprintf('Elapsed time : %6.3g s', etime(clock, starttime)));
else
    % Benchmark mode
    w = floor(params.g_maxiter/5);
    log = zeros(w, 7);
    log(:, [1 4]) = Inf;
    testmaxitersfive = 0;
    bestvalue = Inf;
    best = [];
    alltimes = zeros(params.test_repeat, 1);
    allvalues = zeros(params.test_repeat, 1);
    for iter = 1:params.test_repeat,
        [testbestvalue, testbest, testlog] = genetic(@volumeOBB, @nelderMeadBreed, data, params);
        if (testbestvalue < bestvalue),
            bestvalue = testbestvalue;
            best = testbest;
        end
        nbiters = size(testlog, 1);
        nbitersfive = floor(nbiters/5);
        testmaxitersfive = max(testmaxitersfive, nbitersfive);
        testlog = cell2mat(testlog(5:5:nbiters, 1:2));
        log(1:nbitersfive, 7) = log(1:nbitersfive, 7) + 1;
        log(1:nbitersfive, [3 6]) = log(1:nbitersfive, [3 6]) + testlog;
        log(1:nbitersfive, [2 5]) = max(log(1:nbitersfive, [2 5]), testlog);
        log(1:nbitersfive, [1 4]) = min(log(1:nbitersfive, [1 4]), testlog);
        alltimes(iter) = testlog(end, 1);
        allvalues(iter) = testlog(end, 2);
        elapsed = floor(etime(clock, starttime));
        emin = floor(elapsed/60);
        esec = mod(elapsed, 60);
        ehour = floor(emin/60);
        emin = mod(emin, 60);
        disp(sprintf('Test # %5d/%5d - Finished in %4d iterations with value %17.15g - Cumulative elapsed time : %3dh%2dm%2ds', iter, params.test_repeat, nbiters, testbestvalue, ehour, emin, esec));
    end
    log(:, [3 6]) = log(:, [3 6]) ./ log(:, [7 7]);
    percentiles = [5 10 25 50 75 90 95];
    log = [log(1:testmaxitersfive, :) ; min(alltimes), max(alltimes), mean(alltimes), min(allvalues), max(allvalues), mean(allvalues), std(allvalues) ; prctile(alltimes, percentiles) ; prctile(allvalues, percentiles)];
end

