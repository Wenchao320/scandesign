% script gmj_genetic_test.m
%
% copyright 2017, Mingjie Gao, university of michigan

warning('off');
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

% header options
bool.sv = 0;                                          % save optimized scan design

% construct parameter initialization
rng.de.tr = [17.5 Inf];
rng.sp.tr = [11.8 Inf];
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
  
% gradient function options
subArg.grad = {...
  'x.ff.minmax', [0.03 0.21],...
  'x.ff.nsamp', 5,...
  'x.T1f.nsamp', 1,...
  'x.T1s.nsamp', 1,...
  'x.T2f.nsamp', 1,...
  'x.T2s.nsamp', 1,...
  'nu.kap.nsamp', 3};

% fmincon options
fminconArg = {...
  'boxcon.de.tr', col(rng.de.tr),...
  'boxcon.sp.tr', col(rng.sp.tr),...
  'boxcon.de.aex', col(rng.de.aex),...
  'boxcon.sp.aex', col(rng.sp.aex),...
  'lincon.tr', lincon.tr,...
  'fmincon.tolFun', 1e-7,...
  'fmincon.tolX', 1e-7,...
  'fmincon.alg', 'active-set',...
  'fmincon.disp', 'off',...
  'fmincon.maxIter', 300};

% genetic algorithm
Cde = 3; % will double
Csp = 3;
M = 70;  % population size
Pc = 0.8; % probability of crossover
Pm = 0.2; % probability of mutation
% Pp = 0.3; % probability of copy
G = 80; % total generation

opt.loadpopu = 0;
opt.savepopu = 1;

% initialization of population
if opt.loadpopu==0
    fopt = 100000.0;
    popu = cell(1, M);
    i = 0;
    while i < M
        try
            i = i + 1;
            popu{i}.P.de.tr = rng.de.tr(1) * ones(Cde,1);
            popu{i}.P.sp.tr = rng.sp.tr(1) * ones(Csp,1);
            popu{i}.P.de.aex = col(linspace(randi(60), randi(60), Cde)) * (pi/180);  
            popu{i}.P.sp.aex = col(linspace(randi(40), randi(40), Csp)) * (pi/180);
            popu{i}.f = dess_spgr_2comp_cost(popu{i}.P, subArg.cost{:});
            if isnan(popu{i}.f) || popu{i}.f < 0
                throw(MException('cost value NaN or < 0. ', 'cost value NaN or < 0. '));
            end
            if popu{i}.f < fopt
                fopt = popu{i}.f;
                P.opt = popu{i}.P;
            end
        catch
            i = i - 1;
        end
    end
else
    load(sprintf('popu%dde%dsp%d.mat', M, Cde, Csp), 'popu', 'fopt', 'P');
end
% let's go
g = 0;
tic;
while g < G
    g = g + 1;
    fprintf('\n---------------Generation %d. fopt = %0.6f and rstd = %0.4f.-------------------\nCost functions are:\n', g, fopt, sqrt(fopt) / mean([0.03 0.21]));
    fsum = 0;
    for i = 1:M
        fprintf(' %12.6f', popu{i}.f);
        if popu{i}.f < fopt
            fopt = popu{i}.f;
            P.opt = popu{i}.P;
        end
        fsum  = fsum + popu{i}.f;
        if mod(i,10)==0
            fprintf('\n');
        end
    end
    fprintf('\n');
    
    for i = 1:M
        popu{i}.finv = (fsum - popu{i}.f) / fsum / (M - 1);
        if i == 1
            popu{i}.cummufinv = popu{i}.finv;
        else 
            popu{i}.cummufinv = popu{i - 1}.cummufinv + popu{i}.finv;
        end
    end
    newpopu = cell(1, M);
    
    i = 0;
    while i < M
        try
            i = i + 1;
            % select two individuals
            r = rand;
            for i1 = 1:M
                if r < popu{i1}.cummufinv
                    break;
                end
            end
            if i1 == M 
                i1 = 1;
            else
                i1 = i1 + 1;
            end
            r = rand;
            for i2 = 1:M
                if r < popu{i2}.cummufinv
                    break;
                end
            end
            fprintf('\nGenerate new generation #%d with i1 = %d, i2 = %d. ----> ', i, i1, i2);

            r = rand;
            % crossover
            if r < Pc
                fprintf('crossover');
                deCrossPoint = randi(Cde + 1) - 1;
                spCrossPoint = randi(Csp + 1) - 1;
                newpopu{i}.P.de.aex = [popu{i1}.P.de.aex(1:deCrossPoint); popu{i2}.P.de.aex(deCrossPoint + 1:Cde)];
                newpopu{i}.P.sp.aex = [popu{i1}.P.sp.aex(1:spCrossPoint); popu{i2}.P.sp.aex(spCrossPoint + 1:Csp)];
                newpopu{i}.P.de.tr = [popu{i1}.P.de.tr(1:deCrossPoint); popu{i2}.P.de.tr(deCrossPoint + 1:Cde)];
                newpopu{i}.P.sp.tr = [popu{i1}.P.sp.tr(1:spCrossPoint); popu{i2}.P.sp.tr(spCrossPoint + 1:Csp)];
                if sum(newpopu{i}.P.de.tr(:)) + sum(newpopu{i}.P.sp.tr(:)) > lincon.tr
                    throw(MException('violate time constraint. ', 'violate time constraint. '));
                end
                newpopu{i}.f = dess_spgr_2comp_cost(newpopu{i}.P, subArg.cost{:});
                if isnan(newpopu{i}.f) || newpopu{i}.f < 0
                    throw(MException('cost value NaN or < 0. ', 'cost value NaN or < 0. '));
                end
           
            % mutation
            else % if r >= Pc && r < Pc+Pm
                fprintf('mutation');
                newpopu{i}.P = gmj_Popt_wrapper(popu{i1}.P, subArg, fminconArg{:});
                newpopu{i}.f = dess_spgr_2comp_cost(newpopu{i}.P, subArg.cost{:});
                if isnan(newpopu{i}.f) || newpopu{i}.f < 0
                    throw(MException('cost value NaN or < 0. ', 'cost value NaN or < 0. '));
                end
                  
            % copy
%             else
%                 fprintf('copy');
%                 newpopu{i}.P = popu{i1}.P;
%                 newpopu{i}.f = popu{i1}.f;
%                 if isnan(newpopu{i}.f) || newpopu{i}.f < 0
%                     throw(MException('cost value NaN or < 0. ', 'cost value NaN or < 0. '));
%                 end
            end
        catch ME
            i = i - 1;
            fprintf(' --> Fail!');
        end
    end
    popu = newpopu;
end                                      
t = toc;

% calculate final coefficient of variation w.r.t. mean(ff)
% warning('off', 'MATLAB:nearlySingularMatrix');
rstd.opt  = sqrt(fopt)   ./ mean([0.03 0.21]);

% print output
fprintf('\nTotal execution time is %0.2f minutes.\n', t/60);
fprintf('Optimized profile (%d, %d) yields mean ff rstd = %0.4f with f = %0.6f.\n', Cde, Csp, rstd.opt, fopt);

% save population
if opt.savepopu==1
    save(sprintf('popu%dde%dsp%d.mat', M, Cde, Csp), 'popu', 'fopt', 'P');
end

% save design
if bool.sv
  tmp = sprintf('Popt_TRtot-%.1fms', lincon.tr);
  tmp = strrep(tmp, '.', 'p');
  tmp = strcat(tmp, '.mat');
  save(tmp, 'C0', 'rng', 'P*', 'lincon', '*Arg', 't', 'f', 'rstd*', '-mat');
end
