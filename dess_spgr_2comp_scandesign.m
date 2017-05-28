  function [P, t] = dess_spgr_2comp_scandesign(P0, subArg, varargin)
%|function [P, t] = dess_spgr_2comp_scandesign(P0, subArg, varargin)
%|
%|  optimized scan design for 2-compartment parameter estimation
%|    builds scan profile in a greedy fashion
%|    uses dess_spgr_2comp_Popt(...) to optimize each candidate scan profile
%|    uses dess_spgr_2comp_cost(...) to evaluate cost function at iterates
%|    uses dess_spgr_2comp_costgrad(...) to evaluate cost function gradient
%|
%|  inputs
%|    P0        [1x1 struct]    scan parameters
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|    subArg    [1x1 struct]    subfunction options name-value pairs
%|     .cost    {1 2*noptc}       passed to dess_spgr_2comp_cost(...)
%|     .grad    {1 2*noptg}       passed to dess_spgr_2comp_costgrad(...)
%|
%|  options
%|    boxcon    [1x1 struct]    box constraints, [lower; upper]
%|                                1st arg: data type (sp,de)
%|                                2nd arg: scan parameter
%|     .sp.aex  [2]                 spgr flip angle                     def: [pi/180; pi/2]   rad
%|     .de.aex  [2]                 dess flip angle                     def: [pi/180; pi/2]   rad
%|     .sp.tr   [2]                 spgr repetition time                def: [11.8; Inf]      ms
%|     .de.tr   [2]                 dess repetition time                def: [17.5; Inf]      ms
%|    lincon    [1x1 struct]    linear inequality constraints
%|     .aex     [1]               total aex constraint                  def: Inf              rad
%|     .tr      [1]               total tr constraint                   def: sum(P0.*.tr)     ms
%|    fmincon   [1x1 struct]    fmincon options
%|     .wgrad   [char]            include cost function gradient        def: 'on'
%|     .disp    [char]            console feedback option               def: 'iter-detailed'
%|     .alg     [char]            optimization algorithm                def: 'active-set'
%|     .tolFun  [1]               cost function tolerance               def: 1e-8
%|     .tolX    [1]               iterate tolerance                     def: 1e-10
%|     .maxIter [1]               maximum number of iterations          def: 400
%|
%|  outputs
%|    P         [1x1 struct]    optimized scan design
%|                                1st field: data type (sp,de)
%|                                2nd field: scan parameter, as appropriate
%|     .*.tr    [S.*]               repetition times                                          ms
%|     .*.aex   [S.*]               nominal flip angle of excitation                          rad
%|    t         [1]             scan design run time                                          s
%|
%|  written by: gopal nataraj
%|  copyright 2016, university of michigan
%|
%|  version control
%|    1.1       2016-05-09      original
%|    1.2       2016-08-17      changed format of P*, boxcon

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
fprintf('\nOptimizing (%u,%u) dess/spgr scan profile...\n\n', S.de, S.sp);
[P, cost, flag] = ...
  dess_spgr_2comp_Popt(P0, subArg.cost, subArg.grad, subArg.Popt{:});
if flag>=0
  fprintf('\n...success: initial cost reduced from %0.6f to %0.6f.\n',...
    dess_spgr_2comp_cost(P0, subArg.cost{:}),...
    dess_spgr_2comp_cost(P, subArg.cost{:}));
else
  fprintf('\n...fail: flag = %u.\n', flag);
  return;
end

tic;
while true
  % check to see if we can fit another scan
  if (arg.boxcon.de.tr(1)*S.de + arg.boxcon.sp.tr(1)*S.sp) > arg.lincon.tr
    fprintf('\nExit with (%u,%u) dess/spgr scans:\n', S.de, S.sp);
    fprintf('minimum sum(TR) %0.2f ms exceeds constraint of %0.2f ms.\n\n',...
      arg.boxcon.de.tr(1)*S.de + arg.boxcon.sp.tr(1)*S.sp, arg.lincon.tr);
    break;
  end
  
  % try adding a dess scan
  A.de = 1;
  A.sp = 0;
  P0d = addScan(P, A, arg);
  Snxt.de = S.de+A.de;
  Snxt.sp = S.sp+A.sp;
  subArg.Popt = argBuilder(arg);
  
  fprintf('\n\nTry (%u,%u) dess/spgr scan profile...\n', Snxt.de, Snxt.sp);
  [Pd, costd, flagd] = ...
    dess_spgr_2comp_Popt(P0d, subArg.cost, subArg.grad, subArg.Popt{:});
  if flagd>=0
    fprintf('\n...success:\n');
    fprintf('   cost changed from %0.6f to %0.6f;\n',...
      dess_spgr_2comp_cost(P0d, subArg.cost{:}),...
      dess_spgr_2comp_cost(Pd, subArg.cost{:}));
    fprintf('   total time constraint violation changed from %0.2f ms to %0.2f ms.\n',...
      max(0, sum([P0d.de.tr; P0d.sp.tr]) - arg.lincon.tr),...
      max(0, sum([Pd.de.tr; Pd.sp.tr]) - arg.lincon.tr));
  else
    fprintf('\n...fail: flag = %u.\n', flagd);
  end
  
  % try adding an spgr scan
  A.de = 0;
  A.sp = 1;
  P0s = addScan(P, A, arg);
  Snxt.de = S.de+A.de;
  Snxt.sp = S.sp+A.sp;
  subArg.Popt = argBuilder(arg);
  
  fprintf('\n\nTry (%u,%u) dess/spgr scan profile...\n', Snxt.de, Snxt.sp);
  [Ps, costs, flags] = ...
    dess_spgr_2comp_Popt(P0s, subArg.cost, subArg.grad, subArg.Popt{:});
  if flags>=0
    fprintf('\n...success:\n');
    fprintf('   cost changed from %0.6f to %0.6f;\n',...
      dess_spgr_2comp_cost(P0s, subArg.cost{:}),...
      dess_spgr_2comp_cost(Ps, subArg.cost{:}));
    fprintf('   total time constraint violation changed from %0.2f ms to %0.2f ms.\n',...
      max(0, sum([P0s.de.tr; P0s.sp.tr]) - arg.lincon.tr),...
      max(0, sum([Ps.de.tr; Ps.sp.tr]) - arg.lincon.tr));
  else
    fprintf('\n...fail: flag = %u.\n', flags);
  end
  
  % if neither scan succeeds, report and break
  success = (flagd>=0) | (flags>=0);
  if ~success
    fprintf('\nExit at (%u,%u) dess/spgr profile: next opt steps failed.\n', S.de, S.sp);
    break;
    
  % if neither scan reduces cost, report and break
  elseif cost <= costd && cost <= costs
    fprintf('\nExit at (%u,%u) dess/spgr profile: insufficient reduction in cost.\n', S.de, S.sp);
    break;
    
  % if dess reduces cost more, retain dess addition
  elseif costd <= costs
    fprintf('\nProfile increased from (%u,%u) to (%u,%u) scans.\n', S.de, S.sp, S.de+1, S.sp);
    fprintf('Cost reduced from %0.6f to %0.6f.\n', cost, costd);
    fprintf('Scan time changed from %0.2f ms to %0.2f ms.\n\n',...
      sum([P.de.tr; P.sp.tr]), sum([Pd.de.tr; Pd.sp.tr]));
    P = Pd;
    cost = costd;
    S.de = S.de+1;
    
  % if spgr reduces cost more, retain spgr addition
  else
    fprintf('\nProfile increased from (%u,%u) to (%u,%u) scans.\n', S.de, S.sp, S.de, S.sp+1);
    fprintf('Cost reduced from %0.6f to %0.6f.\n', cost, costs);
    fprintf('Scan time changed from %0.2f ms to %0.2f ms.\n\n',...
      sum([P.de.tr; P.sp.tr]), sum([Ps.de.tr; Ps.sp.tr]));
    P = Ps;
    cost = costs;
    S.sp = S.sp+1;
  end
end

% record total run time (s)
t = toc;
  end

  function P = addScan(Pold, A, arg)
%|function P = addScan(Pold, A, arg)
%|
%|  helper function for appending an old profile with new scans
%|
%|  inputs
%|    Pold      [1x1 struct]    old profile
%|    A         [1x1 struct]    append scan number object
%|     .de      [1]               number of dess scans to append
%|     .sp      [1]               number of spgr scans to append
%|    arg       [1x1 struct]    prescribes box constraints
%|
%|  outputs
%|    P         [1x1 struct]    new profile
%|  
%|  version control
%|    1.1       2016-05-09      original
%|    1.2       2016-08-17      changed format of P*
  
% set new flip angles uniform-randomly within bounds
Pnew.de.aex   = arg.boxcon.de.aex(1) + (arg.boxcon.de.aex(2)-arg.boxcon.de.aex(1)) .* rand(A.de,1);
Pnew.sp.aex   = arg.boxcon.sp.aex(1) + (arg.boxcon.sp.aex(2)-arg.boxcon.sp.aex(1)) .* rand(A.sp,1);

% set new TR values to near-minimum values
Pnew.de.tr    = (arg.boxcon.de.tr(1)+0.001) .* ones(A.de,1);
Pnew.sp.tr    = (arg.boxcon.sp.tr(1)+0.001) .* ones(A.sp,1);

% append to old values
P.de.aex  = [Pold.de.aex; Pnew.de.aex];
P.sp.aex  = [Pold.sp.aex; Pnew.sp.aex];
P.de.tr   = [Pold.de.tr; Pnew.de.tr];
P.sp.tr   = [Pold.sp.tr; Pnew.sp.tr];
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