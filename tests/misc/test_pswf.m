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

% And four more:
val = 1.053221995207094811; % WolframAlpha spheroidalPS(0,0,1,0)
pass(4) = abs(feval(pswf(0,1), 0) - val) < 1e-10;
val = 0.079588252262137702; % WolframAlpha spheroidalPS(100,0,1,0)
pass(5) = abs(feval(pswf(100,1), 0) - val) < 1e-10;
val = -0.005959162025821166; % WolframAlpha spheroidalPS(33,0,33,.9) times -1 !
pass(6) = abs(feval(pswf(33,33), .9) - val) < 1e-10;
val = -0.000000509171971545; % WolframAlpha spheroidalPS(1,0,20,-1) times -1 !
pass(7) = abs(feval(pswf(1,20), -1) - val) < 1e-10;

% Check signs (see DLMNF (30.4.2)):
f = pswf(0:20,100);
L = legpoly(0:20);
pass(8) = all(sign(f(0.001)) == sign(L(0.001)));

end
