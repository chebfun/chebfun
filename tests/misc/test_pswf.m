% Test file for pswf.m 

function pass = test_pswf(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

P = pswf(0:9, 4);

% Check orthogonality:
pass(1) = norm(P'*P - diag(2./(1:2:19))) < 1e-10;

% Check number of roots:
pass(2) = numel(roots(P(:,end))) == 9;

% Check a precomputed value:
P4pi05 = -0.335179227182412994123178747; % WolframAlpha spheroidalPS(4,0,pi,0.5)
pass(3) = abs(feval(pswf(4,pi), 0.5) - P4pi05) < 1e-10;

% And another:
val = 0.175738914563712759; % WolframAlpha spheroidalPS(15,0,sqrt(2),1/3)
pass(4) = abs(feval(pswf(15,sqrt(2)), 1/3) - val) < 1e-10;

% Check signs (see DLMNF (30.4.2)):
f = pswf(0:20,100);
L = legpoly(0:20);
pass(5) = all(sign(f(0.001)) == sign(L(0.001)));

end
