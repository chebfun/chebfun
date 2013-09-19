% Test file for @chebfun/feval.m

function pass = test_feval(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebfun.pref;
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Test empty input.
f = chebfun(@(x) x, pref);
fx = feval(f, []);
pass(1) = isempty(fx);

% Test endpoint evaluation.
pref.chebfun.splitting = 0;
f = chebfun(@(x) erf(x), [-1 1], pref);

lval1 = feval(f, 'left');
lval2 = feval(f, 'start');
lval3 = feval(f, '-');
pass(2) = (abs(lval1 - erf(-1)) < 10*f.epslevel*f.vscale ) ...
    && isequal(lval1, lval2, lval3);

rval1 = feval(f, 'right');
rval2 = feval(f, 'end');
rval3 = feval(f, '+');
pass(3) = (abs(rval1 - erf(1)) < 10*f.epslevel*f.vscale ) ...
    && isequal(rval1, rval2, rval3);

% Test bogus input for endpoint specification.
try
    badval = feval(f, '???');
    pass(4) = false;
catch ME
    pass(4) = strcmp(ME.identifier, 'CHEBFUN:feval:strInput');
end

% Test evaluation with and without left/right limit values.
pref.chebfun.splitting = 1;
f = chebfun(@(x) sign(x), [-1 1], pref);

x = [-0.5 ; 0 ; 0.5];
pass(5) = isequal(feval(f, x), [-1 ; 0 ; 1]);

leftlim1 = feval(f, x, 'left');
leftlim2 = feval(f, x, '-');
pass(6) = isequal(leftlim1, [-1 ; -1 ; 1]) && isequal(leftlim1, leftlim2);

rightlim1 = feval(f, x, 'right');
rightlim2 = feval(f, x, '+');
pass(7) = isequal(rightlim1, [-1 ; 1 ; 1]) && isequal(rightlim1, rightlim2);

% Test bogus input for limit specification string.
try
    badval = feval(f, x, '???');
    pass(8) = false;
catch ME
    pass(8) = strcmp(ME.identifier, 'CHEBFUN:feval:leftrightchar');
end

try
    badval = feval(f, x, []);
    pass(9) = false;
catch ME
    pass(9) = strcmp(ME.identifier, 'CHEBFUN:feval:leftright');
end

% Spot-check a few functions.
pref.chebfun.splitting = 0;
f_exact = @(x) exp(x) - 1;
f = chebfun(f_exact, [-1 1], pref);
x = xr;
pass(10) = (norm(feval(f, x) - f_exact(x), inf) < 10*f.epslevel*f.vscale);

pref.chebfun.splitting = 1;
f_exact = @(x) abs(x)./(1 + x.^2);
f = chebfun(f_exact, [-2 7], pref);
x = 4.5*xr + 2.5;
pass(11) = (norm(feval(f, x) - f_exact(x), inf) < 10*f.epslevel*f.vscale);

pref.chebfun.splitting = 0;
f_exact = @(x) cos(1e4*x);
f = chebfun(f_exact, [1 5], pref);
x = 2*xr + 3;
pass(12) = (norm(feval(f, x) - f_exact(x), inf) < 10*f.epslevel*f.vscale);

pref.chebfun.splitting = 0;
z = exp(2*pi*1i/6);
f_exact = @(t) sinh(t*z);
f = chebfun(f_exact, [-1 1], pref);
x = xr;
pass(13) = (norm(feval(f, x) - f_exact(x), inf) < 10*f.epslevel*f.vscale);

% Check row vector and matrix input.
pref.chebfun.splitting = 0;
f_exact = @(x) cos(x - 0.2);
f = chebfun(f_exact, [-1 1], pref);

x = xr;
err = feval(f, x.') - f_exact(x.');
pass(14) = isequal(size(err), [1 1000]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

x_mtx = reshape(x, [100 10]);
err = feval(f, x_mtx) - f_exact(x_mtx);
pass(15) = isequal(size(err), [100 10]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

x_3mtx = reshape(x, [5 20 10]);
err = feval(f, x_3mtx) - f_exact(x_3mtx);
pass(16) = isequal(size(err), [5 20 10]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

x_4mtx = reshape(x, [5 4 5 10]);
err = feval(f, x_4mtx) - f_exact(x_4mtx);
pass(17) = isequal(size(err), [5 4 5 10]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

% Check behavior for transposed chebfuns.
f.isTransposed = 1; % [TODO]:  Replace with call to transpose().

err = feval(f, x_mtx) - f_exact(x_mtx).';
pass(18) = isequal(size(err), [10 100]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

err = feval(f, x_3mtx) - permute(f_exact(x_3mtx), [2 1 3]);
pass(19) = isequal(size(err), [20 5 10]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

f.isTransposed = 0; % [TODO]:  Replace with call to transpose().

% Check evaluation at points just outside the domain.
x = [(-1 - 1e-6) ; 1e-6 + 1e-6i ; (1 + 1e-6)];
err = feval(f, x.') - f_exact(x.');
pass(20) = norm(err(:), inf) < 10*f.epslevel*f.vscale;

% Check operation for array-valued chebfuns.
pref.chebfun.splitting = 1;
f_exact = @(x) [sin(x) abs(x) exp(1i*x)];
f = chebfun(f_exact, [], pref);

x = xr;
err = feval(f, x) - f_exact(x);
pass(21) = all(max(abs(err)) < 10*f.epslevel*f.vscale);

x_mtx = reshape(x, [100 10]);
err = feval(f, x_mtx) - f_exact(x_mtx);
pass(22) = isequal(size(err), [100 30]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

x_3mtx = reshape(x, [5 20 10]);
err = feval(f, x_3mtx) - f_exact(x_3mtx);
pass(23) = isequal(size(err), [5 60 10]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

x_4mtx = reshape(x, [5 4 5 10]);
err = feval(f, x_4mtx) - f_exact(x_4mtx);
pass(24) = isequal(size(err), [5 12 5 10]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

f.isTransposed = 1;
err = feval(f, x_4mtx) - permute(f_exact(x_4mtx), [2 1 3 4]);
pass(25) = isequal(size(err), [12 5 5 10]) && (norm(err(:), inf) < ...
    10*f.epslevel*f.vscale);

err = feval(f, 'left') - [sin(-1) 1 exp(-1i)].';
pass(26) = norm(err(:), inf) < 10*f.epslevel*f.vscale;

err = feval(f, 'right') - [sin(1) 1 exp(1i)].';
pass(27) = norm(err(:), inf) < 10*f.epslevel*f.vscale;

% Test feval(chebfun, location, 'left'/'right'):
f = chebfun({@(x) [x, 2*x], @(x) [3*x, 4*x]}, [0 1 2]);
pass(28) = norm(feval(f, 1, 'left') - [1, 2], inf) < 10*f.epslevel*f.vscale;
pass(29) = norm(feval(f, 1, 'right') - [3, 4], inf) < 10*f.epslevel*f.vscale;

% Check for dimension mismatch bug when evaluating an array-valued chebfun on a
% vector which contains breakpoint values, repeated several times.
f_exact = @(x) [sin(x) cos(x)];
f = chebfun(f_exact, [-1 1]);
x = [-1 ; 0.5 ; 1];
err = feval(f, x) - f_exact(x);
pass(30) = all(max(abs(err)) < 10*f.epslevel*f.vscale);

end
