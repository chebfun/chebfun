% Test file for @deltafun/cumsum.m

function pass = test_cumsum(pref)

if (nargin < 1)
    pref = chebfunpref();
end

%%
tol = pref.deltaPrefs.deltaTol;

d = deltafun();
F = cumsum(d);
pass(1) = isempty(F);

%%
f = fun.constructor(@(x) exp(x));
mag = rand(5,5);
loc = sort(rand(1,5));

d = deltafun(f, mag, loc);
[F, rval] = cumsum(d);
pass(2) = iscell(F);
%%
f = fun.constructor(@(x) sin(pi*x));
d = deltafun( f, [-1, 1], [-1, 1]);
[F, rJump] = cumsum(d);

pass(3) = ~isa(F, 'deltafun') && feval(F, -1) == -1 ...
    && rJump == 1 && feval(F, 1) == -1; 

%% A test case based on an example by LNT:
n = 6;
x = chebfun('x',[0 n]);
f = 0.5*sin(x);
for j = 1:n-1
  f = f + randn*dirac(x-j);
end
pass(4) = norm(diff(cumsum(f)) - f) < tol;

%%
x = chebfun('x');
f = dirac(x-.5) + dirac(x) + dirac(x+.5) + heaviside(x);
pass(5) = norm(diff(cumsum(f)) - f) < tol;

%% A test case provided by LNT:
x = chebfun('x'); 
f = sign(x)+sign(x-.5); 
f2 = f(-1)+cumsum(diff(f)); 
pass(6) = norm(f-f2) < tol;
end