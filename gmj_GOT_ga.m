  function [P, fval, flag, output, popu] = gmj_GOT_ga(S, costArg, boxconArg, opt)

% default box constraints
arg.boxcon.sp.aex   = col([1 90]) * pi/180;                             % rad
arg.boxcon.de.aex   = col([1 90]) * pi/180;                             % rad
arg.boxcon.sp.tr    = col([11.8 Inf]);                                  % ms
arg.boxcon.de.tr    = col([17.5 Inf]);                                  % ms

% default total linear constraints
arg.lincon.aex = Inf;                                                   % rad
arg.lincon.tr = Inf;                                                    % ms
  
% substitute varargin values as appropriate
arg = vararg_pair(arg, boxconArg);

% check to ensure time constraint is finite and initially feasible
% check whether scan data is enough
if isinf(arg.lincon.tr)
  error('Must explicitly set finite time constraint!');
elseif (arg.boxcon.de.tr(1)*S.de + arg.boxcon.sp.tr(1)*S.sp) > arg.lincon.tr
  error('Time constaint is initially infeasible.');
elseif S.de * 2 + S.sp < 6
    error('Scan data is < 6.');
end

% convert linear constraints
A = [ones(1,S.de+S.sp) zeros(1,S.de+S.sp);...
     zeros(1,S.de+S.sp) ones(1,S.de+S.sp)];                             % [2 2(S.de+S.sp)]
b = [arg.lincon.aex; arg.lincon.tr];                                    % [2 1]

% number of variables
nvars = (S.de + S.sp) * 2;

% no equality constraints
Aeq = [];
beq = [];

% convert box constraints
lb = [...
  ones(S.de,1)*arg.boxcon.de.aex(1);...
  ones(S.sp,1)*arg.boxcon.sp.aex(1);...
  ones(S.de,1)*arg.boxcon.de.tr(1);...
  ones(S.sp,1)*arg.boxcon.sp.tr(1)];
ub = [...
  ones(S.de,1)*arg.boxcon.de.aex(2);...
  ones(S.sp,1)*arg.boxcon.sp.aex(2);...
  ones(S.de,1)*arg.boxcon.de.tr(2);...
  ones(S.sp,1)*arg.boxcon.sp.tr(2)];

% no nonlinear constraints
nonlcon = [];

% anonymous function to handle extra arguments
fn = @(Pv) gmj_GOT_ga_helper(Pv, S, costArg);

% constrained optimization
[Pv, fval, flag, output, popu] = ga(fn, nvars, A, b, Aeq, beq, lb, ub, nonlcon, opt);  % [2(S.de+S.sp) 1]
        
% convert Pv vec -> P struct for output
P.de.aex  = Pv(1:S.de);
P.sp.aex  = Pv(S.de+1:S.de+S.sp);
P.de.tr   = Pv(S.de+S.sp+1:2*S.de+S.sp);
P.sp.tr   = Pv(2*S.de+S.sp+1:end);
end
  
  
  function cost = gmj_GOT_ga_helper(Pv, S, costArg)

Pv = col(Pv);
% convert Pv vec -> P struct
P.de.aex = Pv(1:S.de,1);
P.sp.aex = Pv(S.de+1:S.de+S.sp);
P.de.tr   = Pv(S.de+S.sp+1:2*S.de+S.sp);
P.sp.tr   = Pv(2*S.de+S.sp+1:end);

% call cost function
try
    cost = dess_spgr_2comp_cost(P, costArg{:});
catch ME
    cost = 1e4;
end
end
