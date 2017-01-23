function pass = test_linopAdjoint(pref)

% Get preferences:
if ( nargin < 1 )
    pref = cheboppref();
end
tol = pref.bvpTol; 

%% Building blocks
dom = [-1 1];
I = operatorBlock.eye(dom);
D = operatorBlock.diff(dom);
Z = operatorBlock.zeros(dom);
x = chebfun('x', dom);
c = sin(pi*x.^2);
C = operatorBlock.mult(c);   
z = functionalBlock.zero(dom);
E = functionalBlock.eval(dom);
El = E(dom(1));
Er = E(dom(end));


%% Self-adjoint operator (dirichlet BCs): L*u = -u'' + sin(pi*x.^2).*u
u = (x-dom(1)).*(x-dom(end)).*exp(x); % function that satisfies BCs
nrm = norm(u);
L = linop(-D^2+C);
L = addbc(L,El,0);
L = addbc(L,Er,0);
[Ls,op,bcOpL,bcOpR,bcOpM] = linopAdjoint(L,'bvp');
% test action of operator and function handle
pass(1) = norm(-diff(u,2)+c.*u - Ls*u) < nrm*tol;
pass(2) = norm(-diff(u,2)+c.*u - op(x,u)) < nrm*tol;
% test boundary conditions and functionals
pass(3) = norm(Ls.constraint.functional*u) < nrm*tol;
pass(4) = feval(bcOpL(u),dom(1)) < nrm*tol;
pass(5) = feval(bcOpR(u),dom(end)) < nrm*tol;
pass(6) = isempty(bcOpM);
% check constraint values
pass(7) = norm(Ls.constraint.values) == 0;
% check commutator
pass(8) = abs(u'*(L*u) - (Ls*u)'*u) < nrm*tol;


%% Self-adjoint operator (periodic BCs): L*u = -u'' + sin(pi*x.^2).*u
u = (x-dom(1)).*(x-dom(end)).*exp(x); % function that satisfies BCs
nrm = norm(u);
L = linop(-D^2+C);
[Ls,op,bcOpL,bcOpR,bcOpM] = linopAdjoint(L,'periodic');
% test action of operator and function handle
pass(9) = norm(-diff(u,2)+c.*u - Ls*u) < nrm*tol;
pass(10) = norm(-diff(u,2)+c.*u - op(x,u)) < nrm*tol;
% test boundary conditions and functionals
pass(11) = isempty(Ls.constraint.functional);
pass(12) = isempty(bcOpL);
pass(13) = isempty(bcOpR);
pass(14) = strcmp(bcOpM,'periodic');
% check constraint values
pass(15) = isempty(Ls.constraint.values);
% check commutator
pass(16) = abs(u'*(L*u) - (Ls*u)'*u) < nrm*tol;


%% Operator with no BCs: L*u = u'
u = exp(x); % test function for L
v = (x-dom(1)).*(x-dom(end)).*exp(x); % test function for Ls
nrm = norm(u)+norm(v);
L = linop(D);
[Ls,op,bcOpL,bcOpR,bcOpM] = linopAdjoint(L,'bvp');
% test action of operator and function handle
pass(17) = norm(-diff(v) - Ls*v) < nrm*tol;
pass(18) = norm(-diff(v) - op(x,v)) < nrm*tol;
% test boundary conditions and functionals
pass(19) = norm(Ls.constraint.functional*v) < nrm*tol;
pass(20) = norm(feval(bcOpL(v),dom(1))) < nrm*tol;
pass(21) = norm(feval(bcOpR(v),dom(end))) < nrm*tol;
pass(22) = isempty(bcOpM);
% check constraint values
pass(23) = norm(Ls.constraint.values) == 0;
% check commutator
pass(24) = abs(v'*(L*u) - (Ls*v)'*u) < nrm*tol;


%% Block operator (IVP): L(u1,u2) = [ u1'+u2; u1+u2' ]
u1 = (x-dom(1)).*exp(x); u2 = (x-dom(1)).*sin(x); % test functions for L
v1 = (x-dom(end)).*exp(x); v2 = (x-dom(end)).*sin(x); % test functions for Ls
nrm = norm(u1)+norm(u2)+norm(v1)+norm(v2);
L = linop([D,I;I,D]);
L = addbc(L,[El,z],0);
L = addbc(L,[z,El],0);
[Ls,op,bcOpL,bcOpR,bcOpM] = linopAdjoint(L,'bvp');
% test action of operator and function handle
pass(25) = norm([-diff(v1)+v2;v1-diff(v2)] - Ls*[v1;v2]) < nrm*tol;
pass(26) = norm([-diff(v1)+v2;v1-diff(v2)] - op(x,v1,v2)) < nrm*tol;
% test boundary conditions and functionals
pass(27) = norm(Ls.constraint.functional*[v1;v2]) < nrm*tol;
pass(28) = isempty(bcOpL);
pass(29) = norm(feval(bcOpR(v1,v2),dom(end))) < nrm*tol;
pass(30) = isempty(bcOpM);
% check constraint values
pass(31) = norm(Ls.constraint.values) == 0;
% check commutator
pass(32) = abs([v1;v2]'*(L*[u1;u2]) - (Ls*[v1;v2])'*[u1;u2]) < nrm*tol;


%% Block operator (FVP): L(u1,u2) = [ u1'-u2'; u1+u2' ]
u1 = (x-dom(end)).*exp(x); u2 = (x-dom(end)).*sin(x); % test functions for L
v1 = (x-dom(1)).*exp(x); v2 = (x-dom(1)).*sin(x); % test functions for Ls
nrm = norm(u1)+norm(u2)+norm(v1)+norm(v2);
L = linop([D,-D;I,D]);
L = addbc(L,[Er,z],0);
L = addbc(L,[z,Er],0);
[Ls,op,bcOpL,bcOpR,bcOpM] = linopAdjoint(L,'bvp');
% test action of operator and function handle
pass(33) = norm([-diff(v1)+v2;diff(v1)-diff(v2)] - Ls*[v1;v2]) < nrm*tol;
pass(34) = norm([-diff(v1)+v2;diff(v1)-diff(v2)] - op(x,v1,v2)) < nrm*tol;
% test boundary conditions and functionals
pass(35) = norm(Ls.constraint.functional*[v1;v2]) < nrm*tol;
pass(36) = norm(feval(bcOpL(v1,v2),dom(1))) < nrm*tol;
pass(37) = isempty(bcOpR);
pass(38) = isempty(bcOpM);
% check constraint values
pass(39) = norm(Ls.constraint.values) == 0;
% check commutator
pass(40) = abs([v1;v2]'*(L*[u1;u2]) - (Ls*[v1;v2])'*[u1;u2]) < nrm*tol;

end
