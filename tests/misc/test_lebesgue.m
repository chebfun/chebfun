% Test file for lebesgue.m.
% Based on the test from Chebfun v4.

function pass = test_lebesgue(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

tol = pref.eps;

[L, C] = lebesgue(chebpts(3));
pass(1) = abs(C - 5/4) < 10*tol;

[L, C] = lebesgue(legpts(3), [-1,1]);
pass(2) = abs(C - 7/3) < 11*tol;

[L, C] = lebesgue(linspace(5, 9, 3), 5, 9);
pass(3) = abs(C - 5/4) < 10*tol;

L = lebesgue([1 2], [0, 7]);
pass(4) = abs(norm(L, inf) - 11) < 100*tol;

[L, C] = lebesgue([-1 0 1], [-1, 2]);
pass(5) = norm([7 7] - [C max(L)], inf) < 100*tol;

s = linspace(.25, 1, 17);
s = [-s s];
[L, C] = lebesgue(s);
pass(6) = min(L) > .999;

% This test only applies to chebtechs. See #740.
if ( isa(pref.tech(), 'chebtech') )
    f = lebesgue(linspace(-1, 1, 40));
    pass(7) = length(f) < 1560;
else 
    pass(7) = true;
end


end
