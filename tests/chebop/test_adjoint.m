function pass = test_adjoint(pref)

% Get preferences:
if ( nargin < 1 )
    pref = cheboppref();
end
tol = pref.bvpTol; 

%% Building blocks
dom = [-1 1];
x = chebfun('x', dom);
c = sin(pi*x.^2);
o = 0*x+1;


%% Self-adjoint operator (dirichlet BCs): L*u = -u'' + sin(pi*x.^2).*u
u = (x-dom(1)).*(x-dom(end)).*exp(x); % function that satisfies BCs
nrm = norm(u);
clear L; L = chebop(dom);
L.op = @(x,u) -diff(u,2) + c.*u;
L.lbc = 0;
L.rbc = 0;
Ls = adjoint(L);
% test action of operator and function handle
pass(1) = norm(-diff(u,2)+c.*u-Ls(u)) < nrm*tol;
% test boundary conditions
pass(2) = feval(Ls.lbc(o),dom(1));
pass(3) = feval(Ls.rbc(o),dom(end));
% check commutator
pass(4) = abs(u'*(L*u) - (Ls*u)'*u) < nrm*tol;


%% Self-adjoint operator (periodic BCs): L*u = -u'' + sin(pi*x.^2).*u
u = (x-dom(1)).*(x-dom(end)).*exp(x); % function that satisfies BCs
nrm = norm(u);
clear L; L = chebop(dom);
L.op = @(x,u) -diff(u,2) + c.*u;
L.bc = 'periodic';
Ls = adjoint(L);
% test action of operator and function handle
pass(5) = norm(-diff(u,2)+c.*u - Ls*u) < nrm*tol;
% test boundary conditions
pass(6) = strcmp(Ls.bc,'periodic');
% check commutator
pass(7) = abs(u'*(L*u) - (Ls*u)'*u) < nrm*tol;


%% Operator with no BCs: L*u = u'
u = exp(x); % test function for L
v = (x-dom(1)).*(x-dom(end)).*exp(x); % test function for Ls
nrm = norm(u)+norm(v);
clear L; L = chebop(dom);
L.op = @(x,u) diff(u);
Ls = adjoint(L);
% test action of operator and function handle
pass(8) = norm(-diff(v) - Ls*v) < nrm*tol;
% test boundary conditions 
pass(9) = feval(Ls.lbc(o),dom(1));
pass(10) = feval(Ls.rbc(o),dom(end));
% check commutator
pass(11) = abs(v'*(L*u) - (Ls*v)'*u) < nrm*tol;


%% Block operator (IVP): L(u1,u2) = [ u1'+u2; u1+u2' ]
u1 = (x-dom(1)).*exp(x); u2 = (x-dom(1)).*sin(x); % test functions for L
v1 = (x-dom(end)).*exp(x); v2 = (x-dom(end)).*sin(x); % test functions for Ls
nrm = norm(u1)+norm(u2)+norm(v1)+norm(v2);
clear L; L = chebop(dom);
L.op = @(x,u1,u2) [diff(u1)+u2;u1+diff(u2)];
L.lbc = @(u1,u2) [u1;u2];
Ls = adjoint(L);
% test action of operator and function handle
pass(12) = norm([-diff(v1)+v2;v1-diff(v2)] - Ls*[v1;v2]) < nrm*tol;
% test boundary conditions
pass(13) = min(feval(Ls.rbc(o,o),dom(end)));
% check commutator
pass(14) = abs([v1;v2]'*(L*[u1;u2]) - (Ls*[v1;v2])'*[u1;u2]) < nrm*tol;


%% Block operator (FVP): L(u1,u2) = [ u1'-u2'; u1+u2' ]
u1 = (x-dom(end)).*exp(x); u2 = (x-dom(end)).*sin(x); % test functions for L
v1 = (x-dom(1)).*cos(x); v2 = (x-dom(1)).*sin(x); % test functions for Ls
nrm = norm(u1)+norm(u2)+norm(v1)+norm(v2);
clear L; L = chebop(dom);
L.op = @(x,u1,u2) [diff(u1)-diff(u2);u1+diff(u2)];
L.rbc = @(u1,u2) [u1;u2];
Ls = adjoint(L);
% test action of operator and function handle
pass(15) = norm([-diff(v1)+v2;diff(v1)-diff(v2)] - Ls*[v1;v2]) < nrm*tol;
% test boundary conditions
pass(16) = min(feval(Ls.lbc(o,o),dom(1)));
% check commutator
pass(17) = abs([v1;v2]'*(L*[u1;u2]) - (Ls*[v1;v2])'*[u1;u2]) < nrm*tol;

%% Test L' syntax, and self-adjoint part

Ls1 = adjoint(L);
Ls2 = L';
pass(18) = norm((Ls1*[v1;v2]) - (Ls2*[v1;v2])) < nrm*tol;
pass(19) = abs(imag([v1;v2]'*((L+L')*[v1;v2]))) < nrm*tol;

end
