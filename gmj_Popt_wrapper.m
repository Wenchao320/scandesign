  function [P] = gmj_Popt_wrapper(P0, subArg, varargin)
%|function [P] = gmj_Popt_wrapper(P0, subArg, varargin)
%|
%|  optimized scan design for 2-compartment parameter estimation
%|    builds scan profile in a exhaustive fashion
%|    uses dess_spgr_2comp_Popt(...) to optimize each candidate scan profile
%|    uses dess_spgr_2comp_cost(...) to evaluate cost function at iterates
%|    uses dess_spgr_2comp_costgrad(...) to evaluate cost function gradient

% constant declarations
S.de = length(P0.de.aex);
S.sp = length(P0.sp.aex);

% error checks on TRd, TRs
if length(P0.de.tr) ~= S.de
  error('flipd and TRd of unequal length!');
elseif length(P0.sp.tr) ~= S.sp
  error('flips and TRs of unequal length!');
end

% default initial box constraints
arg.boxcon.sp.aex   = col([1 90]) * pi/180;                             % rad
arg.boxcon.de.aex   = col([1 90]) * pi/180;                             % rad
arg.boxcon.sp.tr    = col([11.8 Inf]);                                  % ms
arg.boxcon.de.tr    = col([17.5 Inf]);                                  % ms

% default total linear constraints
arg.lincon.aex = Inf;                                                   % rad
arg.lincon.tr = sum([P0.sp.tr; P0.de.tr]);                              % ms

% default optimization options
arg.fmincon.wgrad = 'on';     
arg.fmincon.disp = 'iter';
arg.fmincon.alg = 'active-set';
arg.fmincon.tolFun = 1e-8;
arg.fmincon.tolX = 1e-10;
arg.fmincon.maxIter = 400;
  
% substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

% check to ensure time constraint is finite and initially feasible
if isinf(arg.lincon.tr)
  error('Must explicitly set finite time constraint!');
elseif (arg.boxcon.de.tr(1)*S.de + arg.boxcon.sp.tr(1)*S.sp) > arg.lincon.tr
  error('Time constaint is initially infeasible.');
end

% initialization
subArg.Popt = argBuilder(arg);
[P, ~, flag] = ...
  dess_spgr_2comp_Popt(P0, subArg.cost, subArg.grad, subArg.Popt{:});
if flag<0
  fprintf('\n...fail: flag = %u.\n', flag);
  return;
end

  end

  function out = argBuilder(arg)
%|function out = argBuilder(arg)
%|
%|  helper function for constructing an dess_spgr_2comp_Popt(...) options
%|
%|  inputs
%|    S         [1x1 struct]    scan number object
%|     .de      [1]               number of dess scans
%|     .sp      [1]               number of spgr scans
%|    arg       [1x1 struct]    dess_spgr_2comp_scandesign(...) -compatible options
%|
%|  outputs
%|    out       [1x1 struct]    dess_spgr_2comp_Popt(...) -compatible options
%|  
%|  version control
%|    1.1       2016-05-09      original
%|    1.2       2016-08-17      S no longer needed with improved option format

out = {...
  'boxcon.de.aex', arg.boxcon.de.aex,...
  'boxcon.sp.aex', arg.boxcon.sp.aex,...
  'boxcon.de.tr', arg.boxcon.de.tr,...
  'boxcon.sp.tr', arg.boxcon.sp.tr,...
  'lincon.aex', arg.lincon.aex,...
  'lincon.tr', arg.lincon.tr,...
  'fmincon.wgrad', arg.fmincon.wgrad,...
  'fmincon.disp', arg.fmincon.disp,...
  'fmincon.alg', arg.fmincon.alg,...
  'fmincon.tolFun', arg.fmincon.tolFun,...
  'fmincon.tolX', arg.fmincon.tolX,...
  'fmincon.maxIter', arg.fmincon.maxIter};
end