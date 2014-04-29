% Test file for @chebfun/ellipj.m.

function pass = test_ellipj(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points to use as test values.
seedRNG(6178);
xr = 2 * rand(100, 1) - 1;

% Real chebfun u, double m.
pass(1) = do_test_ellipj(@(x) exp(x), 0.6, xr, pref);
pass(2) = do_test_ellipj(@(x) [exp(x) sin(pi*x)], 0.6, xr, pref);

% Real double u, chebfun m.
pass(3) = do_test_ellipj(0.3, @(x) 0.99*abs(x), xr, pref);
pass(4) = do_test_ellipj(0.3, @(x) .05 + abs(.9*[sin(pi*x), cos(pi*x)]), ...
    xr, pref);

% Real chebfun u, chebfun m.
pass(5) = do_test_ellipj(@(x) exp(x), @(x) 0.99*abs(x), xr, pref);
pass(6) = do_test_ellipj(@(x) [exp(x) sin(pi*x)], ...
    @(x) .05 + abs(.9*[sin(pi*x), cos(pi*x)]), xr, pref);

% TODO:  Test the complex case when chebfun.rdivide() is available.

% Check 'tol' option.
u = chebfun(@(x) x, [-1 -0.5 0 0.5 1], pref);
[sn, cn, dn] = ellipj(u, 0.9, 1e-2);
[sn_ex, cn_ex, dn_ex] = ellipj(xr, 0.9, 1e-2);
sn_err = feval(sn, xr) - sn_ex;
cn_err = feval(cn, xr) - cn_ex;
dn_err = feval(dn, xr) - dn_ex;
pass(7) = norm(sn_err(:), inf) < 10*vscale(sn)*epslevel(sn) && ...
          norm(cn_err(:), inf) < 10*vscale(cn)*epslevel(cn) && ...
          norm(dn_err(:), inf) < 10*vscale(dn)*epslevel(dn);
end

% Run one test of chebfun.ellipj().  u_op and m_op are function handles or
% scalars used to supply the u and m inputs to ellipj().  xr is a grid of test
% values in [-1, 1].
function pass = do_test_ellipj(u_op, m_op, xr, pref)
    % Decide if u and m are functions or scalars and convert objects of the
    % former type to chebfuns.
    u_is_fh = isa(u_op, 'function_handle');
    m_is_fh = isa(m_op, 'function_handle');

    if ( u_is_fh )
        u = chebfun(u_op, [-1 -0.5 0 0.5 1], pref);
    else
        u = u_op;
    end

    if ( m_is_fh )
        m = chebfun(m_op, [-1 -0.5 0 0.5 1], pref);
    else
        m = m_op;
    end

    % Call chebfun.ellipj().
    [sn, cn, dn] = ellipj(u, m);

    % Compute the true values.
    if ( u_is_fh && m_is_fh )
        [sn_ex, cn_ex, dn_ex] = ellipj(u_op(xr), m_op(xr));
    elseif ( u_is_fh )
        [sn_ex, cn_ex, dn_ex] = ellipj(u_op(xr), m);
    elseif ( m_is_fh )
        [sn_ex, cn_ex, dn_ex] = ellipj(u, m_op(xr));
    else
        error('Both u and m are scalars.');
    end

    % Compute errors.
    sn_err = feval(sn, xr) - sn_ex;
    cn_err = feval(cn, xr) - cn_ex;
    dn_err = feval(dn, xr) - dn_ex;

    % Compute pass tolerances.
    passtol_sn = 10*vscale(sn)*epslevel(sn);
    passtol_cn = 10*vscale(cn)*epslevel(cn);
    passtol_dn = 10*vscale(dn)*epslevel(dn);

    % Check pass condition.
    pass = norm(sn_err(:), inf) < 10*vscale(sn)*epslevel(sn) && ...
           norm(cn_err(:), inf) < 10*vscale(cn)*epslevel(cn) && ...
           norm(dn_err(:), inf) < 10*vscale(dn)*epslevel(dn);
end
