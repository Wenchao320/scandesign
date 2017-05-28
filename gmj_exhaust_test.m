% script gmj_exhaust_test.m
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
  'fmincon.disp', 'iter',...
  'fmincon.maxIter', 400};

% exhaustive search
fopt = 100000.0;
tic;
for Cde = 0:floor(lincon.tr / rng.de.tr(1))         % max is 108 / 17.5 = 6
    for Csp = 0:floor(lincon.tr / rng.sp.tr(1))     % max is 108 / 11.8 = 9
        if Cde * 2 + Csp < 6
            continue;
        end
        if Cde * rng.de.tr(1) + Csp * rng.sp.tr(1) > lincon.tr
            continue;
        end
        P0.de.tr = rng.de.tr(1) * ones(Cde,1);      % minimum tr
        P0.sp.tr = rng.sp.tr(1) * ones(Csp,1);      % minimum tr
        P0.de.aex = col(linspace(4, 10, Cde)) * (pi/180);   % guess
        P0.sp.aex = flipud(col(linspace(2, 18, Csp))) * (pi/180);   % guess

        % internal optimization
        [P] = gmj_Popt_wrapper(P0, subArg, fminconArg{:});
        f = dess_spgr_2comp_cost(P,  subArg.cost{:});
        fprintf('Scan design (%dDE, %dSP) gives cost function value of %0.6f.\n', Cde, Csp, f);
        if f > 0 && f < fopt
            fopt = f;
            Popt = P;
            C.de = Cde;
            C.sp = Csp;
        end

    end
end                                          
t = toc;

% calculate final coefficient of variation w.r.t. mean(ff)
% warning('off', 'MATLAB:nearlySingularMatrix');
rstd.opt  = sqrt(fopt)   ./ mean([0.03 0.21]);

% print output
fprintf('\nExhaustive search in %0.2f minutes.\n', t/60);
fprintf('Optimized profile (%dDE, %dSP) yields mean ff rstd = %0.4f.\n', C.de, C.sp, rstd.opt);

% save design
if bool.sv
  tmp = sprintf('Popt_%.1fms-%.4f', lincon.tr, rstd.opt);
  tmp = strrep(tmp, '.', 'p');
  tmp = strcat(tmp, '.mat');
  save(tmp);
end
