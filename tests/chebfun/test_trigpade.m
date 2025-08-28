% Test for trigpade.m.
function pass = test_trigpade(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Set an error tolerance.
tol = 1.0e-12;

%% check for emptiness:
p = trigpade(chebfun);
pass(1) = isempty(p);

%%
f = chebfun(@(t) sin(pi*t), 'trig');
m = 1; n = 0;
[p, q, r] = trigpade(f, m, n);
tt = -1+2*rand(100,1);
pass(2) = length(p) <= 2*m+1;
pass(3) = length(q) <= 2*n+1;
pass(4) = norm(p./q - f, inf) < tol;
pass(5) = norm(f(tt) - r(tt), inf) < tol;
[p, q, r, s, t, u, v] = trigpade(f, m, n);
pass(6) = length(p) <= 2*m+1 && length(q) <= 2*n+1;
pass(7) = length(s) <= 2*m+1 && length(t) <= 2*n+1;
pass(8) = length(u) <= 2*m+1 && length(v) <= 2*n+1;
pass(9) = norm(p./q - (s./t + u./v), inf) < tol; 
pass(10) = norm(f(tt) - r(tt), inf) < tol;

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
[p, q] = trigpade(f, m, n);
% approximation has defect, make sure this 
% is detected:
pass(11) = length(p) == 3;
pass(12) = length(q) == 3;
% We should reproduce f:
r = p./q;
pass(13) = norm(f-r, inf) < tol;

%% Check another periodic function 
f = chebfun(@(t) exp(sin(pi*t)), 'trig');
m = 8;
n = 8;
[p, q, r_h] = trigpade(f, m, n);
pass(14) = norm(f-p./q, inf) < tol;
pass(15) = norm(f(tt)-r_h(tt), inf) < tol;


%% Explicitly test the case when m < n
f = chebfun(@(t) 1./(1+sin(pi*(t-.8)/2).^2), 'trig');
m = 0; n = 1;
[p, q] = trigpade(f, m, n);
pass(16) = length(p) == 1;
pass(17) = length(q) == 3;
pass(18) = norm(f-p./q, inf) < tol;

%% Test another function
f = chebfun(@(t) exp(sin(3*pi*t)+cos(pi*t)), 'trig');
m = 5; n = 15;
[p, q] = trigpade(f, m, n);
r = p./q;
pass(19) = compare_coeffs(f, r, m+n, tol);

%% And another function
fh = @(x) (-x+0.7)/(1+0.7).*(x<0.7) + (x-0.7)/(1-0.7).*(x>=0.7);
f = chebfun(fh, 1001, 'trig');
m = 5; n = 5;
[p, q] = trigpade(f, m, n);
r = p./q;
pass(20) = compare_coeffs(f, r, m + n, 100*tol);
pass(21) = length(p) <= 2*m + 1;
pass(22) = length(q) <= 2*n + 1;

%% test m < n again
m = 2; n = 7;
[p, q] = trigpade(f, m, n);
r = p./q;
pass(23) = compare_coeffs(f, r, m + n, 100*tol);
pass(24) = length(p) <= 2*m + 1;
pass(25) = length(q) <= 2*n + 1;


%% Ex 1 p. 381, Baker 1996 Pade Approximant
N = 16;
h = 1/9;k = 1/h+1;
b = 2; c = 1;
c_p = h.^(1:N).';
c_m = k.^(-N:-1).';
c_k = [c_m; b+c; c_p];
f = chebfun(c_k, 'coeffs', 'trig');
m = 1;
n = 2;
[p, q] = trigpade(f, m, n);
r = p./q;
pass(26) = norm(f-r, inf) < tol;
% degree of q should be reduced
pass(27) = length(q) == 2*(n-1)+1; 
pass(28) = length(p) == 2*m+1;


%% Test a complex valued function
f = chebfun(@(t) 1i./(1+25*sin(pi*(t+.5)/2).^2) + exp(sin(3*pi*t)+cos(pi*t)), 'trig');
m = 1; n = 3;
[p, q] = trigpade(f, m, n);
r = p./q;
pass(29) = compare_coeffs(f, r, m+n, tol);
m = 5; n = 2;
[p, q] = trigpade(f, m, n);
r = p./q;
pass(30) = compare_coeffs(f, r, m+n, tol);
end


function out = compare_coeffs(f, g, n, tol)
% Compare 2n+1 coeffs from the middle in
% periodic chebfuns f and g upto a tolerance tol
c1 = f.coeffs;
c2 = g.coeffs;

l1 = (length(c1)-1)/2;
l2 = (length(c2)-1)/2;

% pad with zeros and make them of equal lengths
if ( l1 > l2 )
    c2 = [zeros(l1-l2, 1); c2; zeros(l1-l2, 1)];
else
    c1 = [zeros(l2-l1, 1); c1; zeros(l2-l1, 1)];
end

if ( length(c1) ~= length(c2) )
    error( 'must have same length')
end

% test the coefficients
l = (length(c1)+1)/2;
if (n > l-1)
    out = norm(c1-c2, inf) < tol;
else
    out = norm(c1(l-n:l+n) - c2(l-n:l+n), inf) < tol;
end

end

