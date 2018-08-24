function pass = test_feval( ) 

eps = 1e-10;

% Example 1
F = zeros(10,10,10);
F(2,7,7)=1;
f = ballfun(F);
g = @(r,lam,th)r.*exp(1i*lam).*exp(1i*th);
pass(1) = (feval(f,0.5,1,0.7)==g(0.5,1,0.7));

% Example 2
S = [21,22,23];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th),S);
F = feval(f,[0.5,0.7], 0, pi/2);
exact = [0.5^2, 0.7^2];
pass(2) = max((abs(F(:)-exact(:)))) < eps;

% Example 3
S = [22,23,24];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th),S);
F = feval(f,1, [pi/4,pi/3], pi/2);
exact = [cos(pi/4); cos(pi/3)];
pass(3) = max((abs(F(:)-exact(:)))) < eps;

% Example 3
S = [25,23,20];
f = ballfun(@(r,lam,th)r.^2.*cos(lam).*sin(th),S);
F = feval(f,1, 0, [pi/5, pi/7]);
exact = zeros(1,1,2);
exact(1,1,1) = sin(pi/5);
exact(1,1,2) = sin(pi/7);
pass(3) = max((abs(F(:)-exact(:)))) < eps;

% Example 4
S = [22,23,25];
f = ballfun(@(r,lam,th)r.^3.*sin(lam).*cos(th),S);
r = chebpts(S(1));
lam = pi*trigpts(S(2));
th = pi*trigpts(S(3));
exact = ballfun.coeffs2vals(f.coeffs);
F = feval(f,r,lam,th);
pass(4) = max((abs(F(:)-exact(:)))) < eps;

%% EXTRACT_SPHEREFUN EXAMPLES: 
% Example 1
f = ballfun(@(r,lam,th)cos(lam).*sin(th),S);
g = extract_spherefun(f);
h = spherefun(@(lam,th)cos(lam).*sin(th));
error = coeffs2(g-h,p,n);
pass(1) = max(abs(error(:))) < eps;

% Example 2
f = ballfun(@(r,lam,th)r,[20,21,22]);
g = extract_spherefun(f, 0.5);
h = spherefun(@(lam,th)0.5,p,n);
error = coeffs2(g-h,p,n);
pass(2) = max(abs(error(:))) < eps;

if (nargout > 0)
    pass = all(pass(:));
end
end
