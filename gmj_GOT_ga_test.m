% script gmj_GOT_ga_test.m
%
% copyright 2017, Mingjie Gao, university of michigan
clear all; clc

% setup
if (~exist('irtdir', 'var'))
  curdir = cd('../../irt'); 
  irtdir = pwd;
  setup(); 
  cd(curdir);
end

% add relevant directories
addpath('../model/spgr/');
addpath('../model/dess/');
addpath('../crb/');
addpath('../etc/');

% file output
record = fopen('record.txt', 'w');
fprintf(record, 'This is the record for exhaustive numeration & ga in MATLAB for scan design. \n\n');
fclose(record);

% construct parameter initialization
rng.de.tr = [17.5 500];
rng.sp.tr = [11.8 500];
rng.de.aex = [1 60] * (pi/180);                       % control energy deposition
rng.sp.aex = [1 40] * (pi/180);                       % narrow range to minimize partial spoiling

lincon.tr = (2*3.6 + 4.8) * 9;                      % deoni:11:com time budget at 3.0T = 108

% cost function options
subArg.cost = {...
  'x.ff.minmax', [0.03 0.21],...
  'x.ff.nsamp', 5,...
  'x.T1f.nsamp', 1,...
  'x.T1s.nsamp', 1,...
  'x.T2f.nsamp', 1,...
  'x.T2s.nsamp', 1,...
  'x.kfs.nsamp', 1,...
  'nu.kap.nsamp', 3};

% boxcon options
boxconArg = {...
  'boxcon.de.tr', col(rng.de.tr),...
  'boxcon.sp.tr', col(rng.sp.tr),...
  'boxcon.de.aex', col(rng.de.aex),...
  'boxcon.sp.aex', col(rng.sp.aex),...
  'lincon.tr', lincon.tr};

% hybrid function patternsearch options
patternopt = optimoptions('patternsearch',...
    'FunctionTolerance', 1e-7,...
    'Display', 'iter',...
    'UseParallel', true,...
    'UseCompletePoll', true);

% ga optimization options
opt = optimoptions('ga',...
    'CreationFcn', @gacreationlinearfeasible,...
    'CrossoverFcn', @crossoverintermediate,...
    'MutationFcn', @mutationadaptfeasible,...
    'FunctionTolerance', 1e-7,...
    'MaxGenerations', 500,...
    'MaxStallGenerations', 40,...
    'UseParallel', true, ...
    'UseVectorized', false,...
    'HybridFcn', {@patternsearch,patternopt},...
    'Display', 'iter');

warning('off');

% exhaustive numeration & genetic algorithm
maxCde = floor(lincon.tr / rng.de.tr(1));
maxCsp = floor(lincon.tr / rng.sp.tr(1));
for i1 = 1:maxCde
    for i2 = 1:maxCsp
        C.de = i1;
        C.sp = i2;
        
        if (rng.de.tr(1)*C.de + rng.sp.tr(1)*C.sp) > lincon.tr
          	continue;
        elseif C.de * 2 + C.sp < 6
            continue;
        end
        
        fprintf('With (%dDE, %dSP), let''s call ga. \n', C.de, C.sp);
        tic;
        [P, fval, flag, output, popu] = gmj_GOT_ga(C, subArg.cost, boxconArg, opt);                              
        t = toc;

        % calculate final coefficient of variation w.r.t. mean(ff)
        rstd.opt  = sqrt(fval)   ./ mean([0.03 0.21]);

        record = fopen('record.txt', 'a');
        % print output
        fprintf(record, '\nOptimized profile (%dDE, %dSP) yields mean ff rstd = %0.4f with f = %0.6f. \n',...
            C.de, C.sp, rstd.opt, fval);
        fprintf(record, 'Total execution time is %0.2f minutes.\n', t/60);
        fclose(record);

        % save design
        tmp = sprintf('ga_%dde%dsp-%.1f-%.4f', C.de, C.sp, lincon.tr, rstd.opt);
        tmp = strrep(tmp, '.', 'p');
        tmp = strcat(tmp, '.mat');
        tmp = strcat('./GOT_ga/', tmp);
        save(tmp);
    end
end