% Test file for trigtech/diff.m

function pass = test_diff(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
x = 2 * rand(100, 1) - 1;

testclass = trigtech();

%%
% Spot-check derivatives for a couple of functions.

f = testclass.make({[],[-.25;.75]});
df = diff(f);
df_coeffs_exact = [.25*pi*1i;0];
pass(1) = (norm(df.coeffs-df_coeffs_exact, inf) < eps*norm(df_coeffs_exact,inf));

f = testclass.make(@(x) exp(cos(pi*x)), [], pref);
df = diff(f);
df_exact = @(x) -pi*sin(pi*x).*exp(cos(pi*x));
err = df_exact(x) - feval(df, x);
pass(2) = (norm(err, inf) < 1e3*vscale(df)*eps);
    

a = 10; b = 20;
f = testclass.make(@(x) cos(a*pi*sin(b*pi*x)), [], pref);
df = diff(f);
df_exact = @(x) -pi^2*a*b*cos(b*pi*x).*sin(a*pi*sin(b*pi*x));
err = df_exact(x) - feval(df, x);
pass(3) = (norm(err, inf) < 1e4*vscale(df)*eps);
    

f = testclass.make(@(x) exp(-50*x.^2), [], pref);
df = diff(f);
df_exact = @(x) -100*x.*exp(-50*x.^2);
err = df_exact(x) - feval(df, x);
pass(4) = (norm(err, inf) < 100*vscale(df)*eps);
    

a1 = 4; b1 = 3; a2 = 6; b2 = 4;
f = testclass.make(@(x) cos(a1*pi*sin(b1*pi*x)) + 1i*cos(a2*pi*sin(b2*pi*x)), [], pref);
df = diff(f);
df_exact = @(x) -pi^2*a1*b1*cos(b1*pi*x).*sin(a1*pi*sin(b1*pi*x)) - 1i*pi^2*a2*b2*cos(b2*pi*x).*sin(a2*pi*sin(b2*pi*x));
err = df_exact(x) - feval(df, x);
pass(5) = (norm(err, inf) < 1e3*vscale(df)*eps);
    

%%
% Verify that calling diff() gives the same answer as direct construction.

f = testclass.make(@(x) 1/21/pi*cos(21*pi*x), [], pref);
df = testclass.make(@(x) -sin(21*pi*x), [], pref);
err = diff(f) - df;
pass(6) = (norm(err.coeffs, inf) < 100*vscale(df)*eps);

%%
% Verify basic differentiation rules.

f = testclass.make(@(x) exp(1)-exp(cos(3*pi*x)), [], pref);
df = diff(f);
g = testclass.make(@(x) exp(-sin(2*pi*x)), [], pref);
dg = diff(g);
tol_f = 10*vscale(df)*eps;
tol_g = 10*vscale(dg)*eps;

errfn = diff(f + g) - (df + dg);
err = feval(errfn, x);
pass(7) = (norm(err, inf) < 10*max(tol_f, tol_g));
    
    
errfn = diff(f.*g) - (f.*dg + g.*df);
err = feval(errfn, x);
pass(8) = (norm(err, inf) < length(f)*max(tol_f, tol_g));

const = testclass.make(@(x) ones(size(x)), [], pref);
dconst = diff(const);
err = feval(dconst, x);
pass(9) = (norm(err, inf) == 0);

%%
% Check higher-order derivatives.  (NB:  We relax the tolerance by n + 1
% factors of 10, where n is the number of derivatives taken.)

f = testclass.make(@(x) exp(cos(4*pi*x))-1, [], pref);
df2 = diff(f, 2);
df2_exact = @(x) -16*pi^2*exp(cos(4*pi*x)).*(cos(4*pi*x) + cos(4*pi*x).^2 - 1);
err = df2_exact(x) - feval(df2, x);
pass(10) = (norm(err, inf) < 1e3*vscale(df2)*eps);
    

f = testclass.make(@(x) sin(pi*x), [], pref);
df6 = diff(f, 6);
df6_exact = -pi^6*f;
err = feval(df6_exact,x) - feval(df6, x);
pass(11) = (norm(err, inf) < 100*vscale(df6)*eps);

f = testclass.make(@(x) (1/10/pi)*cos(10*pi*sin(pi*x)), [], pref);
df5 = diff(f, 5);
err = feval(df5,[-1;1]);  % Odd derivatives of this function vanish at +-1
pass(12) = (norm(err, inf) < 1e3*vscale(df5)*eps);
    

%%
% Check operation for array-valued chebtech objects.
f = testclass.make(@(x) [exp(-50*x.^2) sin(4*pi*(x-0.2)) 1i*exp(cos(pi*x))], [], pref);
df = diff(f);
df_exact = @(x) [-100*x.*exp(-50*x.^2) 4*pi*cos(4*pi*(x-0.2)) -pi*1i*sin(pi*x).*exp(cos(pi*x))];
err = feval(df, x) - df_exact(x);
pass(13) = (norm(err(:), inf) < 100*max(vscale(df)*eps));
    

% DIM option.
dim2df = diff(f, 1, 2);
g = @(x) [(sin(4*pi*(x-0.2))-exp(-50*x.^2))  (1i*exp(cos(pi*x))-sin(4*pi*(x-0.2)))];
err = feval(dim2df, x) - g(x);
pass(14) = isequal(size(vscale(dim2df)), [1 2]) && ...
    (norm(err(:), inf) < 100*max(vscale(dim2df)*eps));
    


dim2df2 = diff(f, 2, 2);
g = @(x) exp(-50*x.^2) - 2*sin(4*pi*(x-0.2)) + 1i*exp(cos(pi*x));
err = feval(dim2df2, x) - g(x);
pass(15) = isequal(size(vscale(dim2df2)), [1 1]) && ...
    (norm(err(:), inf) < 1e3*max(vscale(dim2df2)*eps));
    

% DIM option should return an empty trigtech for non-array-valued input.
f = testclass.make(@(x) sin(pi*x));
dim2df = diff(f, 1, 2);
pass(16) = (isempty(dim2df.coeffs));

% even example with complex coefficients
f = testclass.make({[],[1+1i;1-1i]});
df = diff(f);
df_coeffs_exact = [-pi*1i*(1+1i);0;];
pass(17) = (norm(df.coeffs-df_coeffs_exact, inf) < eps*norm(df_coeffs_exact,inf));

end
