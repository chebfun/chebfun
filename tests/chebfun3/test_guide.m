function pass = test_guide(pref)
% Test Chebfun3 guide commands.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e3*pref.cheb3Prefs.chebfun3eps;

% commands from guide1: 
f = chebfun3(@(x,y,z) 1./(1+x.^2+y.^2+z.^2));
pass(1) = abs(f(0, 0.5, 0.5) - 2/3) < tol;

pass(2) = abs(sum3(f) - 4.28685406230184188268) < tol;

pass(3) = abs(mean3(f) - (1/8)*4.28685406230184188268) < tol;

pass(4) = abs(max3(f) - 1) < tol;

f1 = chebfun(@(x) exp(x));
len1D = length(f1);
f3 = chebfun3(@(x,y,z) exp(x));
[m3, n3, p3] = length(f3);
pass(5) = m3 < 2* len1D && n3 == 1 && p3 == 1;

f3 = chebfun3(@(x,y,z) exp(y));
[m3, n3, p3] = length(f3);
pass(6) = m3 == 1 && n3 < 2* len1D && p3 == 1;

f3 = chebfun3(@(x,y,z) exp(z));
[m3, n3, p3] = length(f3);
pass(7) = m3 == 1 && n3 == 1 && p3 < 2* len1D; 

% Is || f * g || <= || f || * || g || ?
f = chebfun3(@(x,y,z) sin(x+y.*z));
g = chebfun3(@(x,y,z) cos(15*exp(z))./(5+x.^3+2*y.^2+z));
pass(8) = max3(f.*g) <= max3(f) * max3(g);

% Test line integration:
curve = chebfun(@(t) [cos(t) sin(t) t/(8*pi)], [0, 8*pi]);
f = chebfun3(@(x,y,z) x+y.*z);
I = integral(f, curve);
exact = -sqrt(1+(8*pi)^2)/(8*pi);
pass(9) = abs(I - exact) < tol;

% Check if the 'trig' flag needs less coefficients for a triply periodic
% function:
ff = @(x,y,z) tanh(3*sin(x))-(sin(y+1/2)).^2+cos(6*z);
dom = [-pi pi -pi pi -pi pi];
fTrig = chebfun3(ff, dom, 'trig');
[m, n, p] = length(fTrig);
fCheb = chebfun3(ff, dom);
[m_fCheb, n_fCheb, p_fCheb] = length(fCheb);
pass(10) = m <= m_fCheb && n <= n_fCheb && p <= p_fCheb;

%% Test Higher-order SVD:
f = chebfun3(@(x,y,z) sin(x+2*y+3*z));
[sv, Score, Scols, Srows, Stubes] = hosvd(f);
sv1 = sv{1};
sv2 = sv{2};
sv3 = sv{3};
% Test decay of singular values:
pass(11) = sv1(2) <= sv1(2);

pass(12) = sv2(2) <= sv2(2);

pass(13) = sv3(2) <= sv3(2);

% Test orthogonality in columns, rows and tubes:
pass(14) = norm(eye(size(Scols,2)) - Scols'*Scols) < tol;

pass(15) = norm(eye(size(Srows,2)) - Srows'*Srows) < tol;

pass(16) = norm(eye(size(Stubes,2)) - Stubes'*Stubes) < tol;

% Test the _all orthogonality_ property of the core tensor:
pass(17) = norm(squeeze(Score(1,:,:) .* Score(2,:,:))) < tol;

pass(18) = norm(squeeze(Score(:,1,:) .* Score(:,2,:))) < tol;

pass(19) = norm(Score(:,:,1) .* Score(:,:,2)) < tol;

end