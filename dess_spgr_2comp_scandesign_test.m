% script dess_spgr_2comp_scandesign_test.m
% header file to test dess_spgr_2comp_scandesign(...)
%
% copyright 2016, gopal nataraj, university of michigan
%
% version control
%   1.1     2016-05-09      original
%   1.2     2016-08-16      changed format of P and grad

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
warning('off');

% construct parameter initialization
C0.de = 3;  % will double                                          % max is 263.7 / 17.5 = 15
C0.sp = 3;                                            % max is 263.7 / 11.8 = 22
rng.de.tr = [17.5 Inf];
rng.sp.tr = [11.8 Inf];
rng.de.aex = [1 60] * (pi/180);                       % control energy deposition
rng.sp.aex = [1 40] * (pi/180);                       % narrow range to minimize partial spoiling

P0.de.aex = col(linspace(5, 20, C0.de)) * (pi/180);
P0.sp.aex = flipud(col(linspace(2, 18, C0.sp)) * (pi/180));
P0.de.tr = rng.de.tr(1) * ones(C0.de,1);              % minimum tr
P0.sp.tr = rng.sp.tr(1) * ones(C0.sp,1);              % minimum tr
% load('Popt_108-0p1230.mat', 'P');
% P0.de.aex = col([P.de.aex; P.de.aex]);
% P0.sp.aex = col([P.sp.aex; P.sp.aex]);
% P0.de.tr = col([P.de.tr; P.de.tr]);
% P0.sp.tr = col([P.sp.tr; P.sp.tr]);

% linear constraints
% lincon.tr = (2*3.5 + 5.4) * 9;                        % deoni:11:com time budget at 1.5T = 111.6
lincon.tr = (2*3.6 + 4.8) * 9;                      % deoni:11:com time budget at 3.0T = 108
% lincon.tr = (17.5 + 11.8) * 9;                        % allows 9 dess,spgr scans at min tr

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
  'fmincon.alg', 'active-set',...
  'fmincon.disp', 'iter',...
  'fmincon.tolFun', 1e-7,...
  'fmincon.tolX', 1e-7,...
  'fmincon.maxIter', 500};

% header options
bool.sv = 0;                                          % save optimized scan design

% scan design
[P, t] = dess_spgr_2comp_scandesign(P0, subArg, fminconArg{:});

% compare initial and final coefficient of variation w.r.t. mean(ff)
% warning('off', 'MATLAB:nearlySingularMatrix');
f.init    = dess_spgr_2comp_cost(P0, subArg.cost{:});
f.opt     = dess_spgr_2comp_cost(P,  subArg.cost{:});
rstd.init = sqrt(f.init)  ./ mean([0.03 0.21]);
rstd.opt  = sqrt(f.opt)   ./ mean([0.03 0.21]);

% print output
fprintf('\nScan designed in %0.2f minutes.\n', t/60);
fprintf('Optimized profile yields mean ff rstd = %0.4f.\n', rstd.opt);

% save design
if bool.sv
  tmp = sprintf('Popt_TRtot-%.1fms', lincon.tr);
  tmp = strrep(tmp, '.', 'p');
  tmp = strcat(tmp, '.mat');
  save(tmp, 'C0', 'rng', 'P*', 'lincon', '*Arg', 't', 'f', 'rstd*', '-mat');
end
