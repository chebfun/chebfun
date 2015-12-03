function pass = test_kronOp(pref) 
% Test the chebfun/kron() command, resulting in OPERATORBLOCK objects. 

if ( nargin < 1 ) 
    pref = chebfunpref(); 
end 

tol = 1000*pref.eps; 
d = [0,1];
x = chebfun(@(x) x, d); 
xx1 = chebpts(200, d, 1);
xx2 = chebpts(200, d, 2);
xxe = trigpts(200, d);
%% Simple outer product, only single column CHEBFUNs involved
f = sin(x);
g = tanh(x);
h = cos(x);
A = kron(f, g', 'op');
Ah = f * (g'*h);
AC = chebmatrix(A);
% Operational form
pass(1) = norm(Ah - A*h) < tol;
% Discrete form
options = cheboppref();
options.discretization = @chebcolloc1;
pass(2) = norm(Ah(xx1) - matrix(AC, 200, options)*h(xx1)) < tol;
options.discretization = @chebcolloc2;
pass(3) = norm(Ah(xx2) - matrix(AC, 200, options)*h(xx2)) < tol;
%% Single column CHEBFUNs, periodic to test trigcolloc as well
f = exp(sin(4*pi*x));
g = tanh(.5*cos(2*pi*x));
h = cos(2*pi*x);
A = kron(f, g', 'op');
Ah = f * (g'*h);
AC = chebmatrix(A);
% Operational form
pass(4) = norm(Ah - A*h) < tol;
% Discrete form
options.discretization = @chebcolloc1;
pass(5) = norm(Ah(xx1) - matrix(AC, 200, options)*h(xx1)) < tol;
options.discretization = @chebcolloc2;
pass(6) = norm(Ah(xx2) - matrix(AC, 200, options)*h(xx2)) < tol;
options.discretization = @trigcolloc;
pass(7) = norm(Ah(xxe) - matrix(AC, 200, options)*h(xxe)) < tol;
%% Multiple columns/rows (taken from v4 test)
f = [ exp(x), tanh(x) ];
g = [ exp(x), x./(1+x.^2) ];
u = x;

A = kron(f,g','op');
AC = chebmatrix(A);
Au = (exp(x) + (1-pi/4)*tanh(x));

% Operational form
pass(8) = norm( Au - A*u ) < tol;
% Discrete form
options.discretization = @chebcolloc1;
pass(9)  = norm( Au(xx1) - matrix(AC, 200, options)*u(xx1) ) < tol;
options.discretization = @chebcolloc2;
pass(10) = norm( Au(xx2) - matrix(AC, 200, options)*u(xx2) ) < tol;

%% Invalid calls to kron with op output:

try
    kron(f, g, 'invalid')
    pass(11) = 0;
catch ME
    pass(11) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:kron:sizes');
end

try
    kron(f.', g, 'op')
    pass(12) = 0;
catch ME
    pass(12) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:kron:columnsAndRows');
end

try
    kron(quasimatrix(f), g', 'op')
    pass(13) = 0;
catch ME
    pass(13) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:kron:quasimatrix');
end

end
