% Test for trigpade.m.
function pass = test_trigpade(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Set an error tolerance.
tol = 1.0e-13;

%% check for emptiness:
p = trigpade(chebfun);
pass(1) = isempty(p);

%%
f = chebfun(@(t) sin(pi*t), 'trig');
[p, q, r] = trigpade(f, 1, 0);
pass(2) = length(p) == 3;
pass(3) = length(q) == 1;
pass(4) = norm(p./q - f, inf) < tol;
[p, q, r, s, t, u, v] = trigpade(f, 1, 0);
pass(5) = length(p) == 3 && length(q) == 1;
pass(6) = length(s) == 3 && length(t) == 1;
pass(7) = length(u) == 3 && length(v) == 1;
pass(8) = norm(p./q - (s./t + u./v), inf) < tol; 

%%
% Ex 1 p. 381, Baker 1996 Pade Approximant
h = 1/9;
k = 1/h;
b = 2;
c = 1;
N = 16;
c_p = h.^(1:N).';
c_m = k.^(-N:-1).';
coeffs = [c_m; b+c; c_p];
f = chebfun(coeffs, 'coeffs', 'trig');
m = 1;
n = 2;
[p, q, rh] = trigpade(f, m, n);
% approximation has defect, make sure this 
% is detected:
pass(9) = length(p) == 3;
pass(10) = length(q) == 3;
% We should reproduce f:
r = p./q;
pass(11) = norm(f-r, inf) < tol;

%% Check another periodic function 
f = chebfun(@(t) exp(sin(pi*t)), 'trig');
m = 15;
n = 3;
[p, q, r_h] = trigpade(f, m, n);
pass(12) = norm(f-p./q, inf) < tol;

end