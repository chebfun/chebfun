% Test file for besselroots.m.

function pass = test_besselroots(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Check residual with besselj.
res = @(v) norm(besselj(v, besselroots(v, 30)), Inf);

pass(1) = res(0) < 1e-12;
pass(2) = res(-1) < 1e-10;
pass(3) = res(2.5) < 1e-10;
pass(4) = res(4) < 1e-10;

% We expect less accuracy in the first few roots in these cases.
res = besselj(6, besselroots(6, 30));
pass(5) = (norm(res(1:6), Inf) < 1e-2) && (norm(res(7:end), Inf) < 1e-8);
res = besselj(7, besselroots(7, 30));
pass(6) = (norm(res(1:6), Inf) < 1e-2) && (norm(res(7:end), Inf) < 1e-8);

end
