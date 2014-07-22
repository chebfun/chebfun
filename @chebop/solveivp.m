function y = solveivp(N, varargin)

disp('Initial value problem detected.')

% Convert to first order format
anonFun = treeVar.toFirstOrder(N.op, N.domain);

%% Obtain information about the initial conditions.
% Begin by evaluation N.lbc with a zero chebfun to pick up the desired values:
cheb0 = chebfun(@(x) 0*x, N.domain);
Nl0 = N.lbc(cheb0);
leftEnd = N.domain(1);
initVals = -[Nl0{1}(leftEnd), Nl0{2}(leftEnd)];

% Create an ODESET struct for specifying tolerance:
opts = odeset('absTol',1e-12,'relTol',1e-12);

[t, y]=chebfun.ode113(anonFun, N.domain, initVals, opts);
% To fit in with CHEBOP backslash, just return the solution (first column of the
% output), rather than all columns (the second column corresponds to the
% first derivative etc.).
y = y(:,1);
end