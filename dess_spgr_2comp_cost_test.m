% script dess_spgr_2comp_cost_test.m
% header file to test dess_spgr_2comp_cost(...)
%
% copyright 2016, gopal nataraj, university of michigan
%
% version control
%   1.1     2016-04-07      original
%   1.2     2016-08-16      changed format of P

% setup
if ~exist('irtdir', 'var')
    curdir = cd('../../irt'); 
    irtdir = pwd;
    setup(); 
    cd(curdir);
end

addpath('../model/spgr/');
addpath('../model/dess/');
addpath('../crb/');
addpath('../etc/');

% construct parameter object
S.de = 3;
S.sp = 3;
load('negvalue-setup.mat','P');
% P.de.aex  = (pi/180) * linspace(1, 180, S.de)';   % rad
% P.sp.aex  = (pi/180) * linspace(1, 180, S.sp)';   % rad
% P.de.tr   = 20 * ones(S.de,1);                    % ms
% P.sp.tr   = 15 * ones(S.sp,1);                    % ms

% evaluate run time for point distribution
tic;
c.test = dess_spgr_2comp_cost(P,...
    'x.ff.nsamp', 5,...
    'x.T1f.nsamp', 1,...
    'x.T1s.nsamp', 1,...
    'x.T2f.nsamp', 1,... 
    'x.T2s.nsamp', 1,...
    'x.kfs.nsamp', 1,...
    'nu.kap.nsamp', 3);
t.test = toc;
fprintf('Small-test: run time = %0.3f ms; expected cost = %0.9f.\n',...
  t.test*1000, c.test);

% evaluate run time for default distribution
tic;
c.def = dess_spgr_2comp_cost(P);
t.def = toc;
fprintf('Default: run time = %0.3f s; expected cost = %0.9f.\n',...
  t.def, c.def);
