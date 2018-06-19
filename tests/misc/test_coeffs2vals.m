function pass = test_coeffs2vals(pref)
% test various coeffs2vals, vals2coeffs, vals2vals, and coeffs2coeffs codes.
if ( nargin == 0 )
    pref = chebfunpref();
end

tol = 1000*eps;

f = chebfun(@exp, pref);
N = length(f);
c_leg = legcoeffs(f);
c_cheb = chebcoeffs(f);
v_leg = f(legpts(N));
v_cheb1 = f(chebpts(N,1));
v_cheb2 = f(chebpts(N,2));

pass(1) = norm(c_leg - chebcoeffs2legcoeffs(c_cheb), inf) < tol;
pass(2) = norm(c_leg - legvals2legcoeffs(v_leg), inf) < tol;
pass(3) = norm(c_leg - chebvals2legcoeffs(v_cheb1, 1), inf) < tol;
pass(4) = norm(c_leg - chebvals2legcoeffs(v_cheb2, 2), inf) < tol;

pass(5) = norm(c_cheb - legcoeffs2chebcoeffs(c_leg), inf) < tol;
pass(6) = norm(c_cheb - legvals2chebcoeffs(v_leg), inf) < tol;
pass(7) = norm(c_cheb - chebvals2chebcoeffs(v_cheb1, 1), inf) < tol;
pass(8) = norm(c_cheb - chebvals2chebcoeffs(v_cheb2, 2), inf) < tol;

pass(9) = norm(v_leg - legcoeffs2legvals(c_leg), inf) < tol;
pass(10) = norm(v_leg - chebcoeffs2legvals(c_cheb), inf) < tol;
pass(11) = norm(v_leg - chebvals2legvals(v_cheb1, 1), inf) < tol;
pass(12) = norm(v_leg - chebvals2legvals(v_cheb2, 2), inf) < tol;

pass(13) = norm(v_cheb1 - legcoeffs2chebvals(c_leg, 1), inf) < tol;
pass(14) = norm(v_cheb1 - chebcoeffs2chebvals(c_cheb,1), inf) < tol;
pass(15) = norm(v_cheb1 - legvals2chebvals(v_leg, 1), inf) < tol;
pass(16) = norm(v_cheb1 - chebvals2chebvals(v_cheb2, 2, 1), inf) < tol;

pass(17) = norm(v_cheb2 - legcoeffs2chebvals(c_leg, 2), inf) < tol;
pass(18) = norm(v_cheb2 - chebcoeffs2chebvals(c_cheb,2), inf) < tol;
pass(19) = norm(v_cheb2 - legvals2chebvals(v_leg, 2), inf) < tol;
pass(20) = norm(v_cheb2 - chebvals2chebvals(v_cheb1, 1, 2), inf) < tol;

end