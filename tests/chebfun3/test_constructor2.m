function pass = test_constructor2( pref ) 
% This tests the chebfun3 constructor. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb3Prefs.chebfun3eps;

% Test building Chebfun3 objects from sample data: 
exactCoeffs = rand(3,3,3);
constructorCoeffs = chebcoeffs3(chebfun3(exactCoeffs, 'coeffs'));
pass(1) = norm(exactCoeffs(:) - constructorCoeffs(:)) < 10*tol;

exactCoeffs = rand(4,4,4);
constructorCoeffs = chebcoeffs3(chebfun3(exactCoeffs, 'coeffs'));
pass(2) = norm(exactCoeffs(:) - constructorCoeffs(:)) < 10*tol;

% Make a chebfun3 based on TRIGTECH by calling the 'trig' flag:
f1 = chebfun3(@(x,y,z) cos(pi*x).*sin(pi*y).*cos(pi*z),'trig'); 
f2 = chebfun3(@(x,y,z) cos(pi*x).*sin(pi*y).*cos(pi*z), ...
    [-1 1 -1 1 -1 1], 'trig');
pass(3) = norm( f1 - f2 ) < tol;

f1 = chebfun3(@(x,y,z) cos(pi*cos(pi*x) + pi*sin(pi*y) + ...
    pi*cos(pi*z)),'trig');
f2 = chebfun3(@(x,y,z) cos(pi*cos(pi*x) + pi*sin(pi*y) + ...
    pi*cos(pi*z)), [-1 1 -1 1 -1 1], 'trig');
pass(4) = ( norm(f1 -f2) < 10*tol );

% Check underlying tech is a TRIGTECH: 
colTech = get(f1.cols.funs{1}, 'tech');
rowTech = get(f1.rows.funs{1}, 'tech');
tubeTech = get(f1.tubes.funs{1}, 'tech');
pass(5) = isa(colTech(), 'trigtech'); 
pass(6) = isa(rowTech(), 'trigtech'); 
pass(7) = isa(tubeTech(), 'trigtech'); 

% Make sure the 'periodic' flag works as well:
f1 = chebfun3(@(x,y,z) cos(pi*cos(pi*x) + pi*sin(pi*y) + ...
    pi*cos(pi*z)),'periodic');
f2 = chebfun3(@(x,y,z) cos(pi*cos(pi*x) + pi*sin(pi*y) + ...
    pi*cos(pi*z)), [-1 1 -1 1 -1 1], 'periodic');
pass(8) = ( norm(f1 -f2) < 10*tol );

% Check underlying tech is a TRIGTECH:
colTech = get(f1.cols.funs{1}, 'tech');
rowTech = get(f1.rows.funs{1}, 'tech');
tubeTech = get(f1.tubes.funs{1}, 'tech');
pass(9) = isa(colTech(), 'trigtech'); 
pass(10) = isa(rowTech(), 'trigtech'); 
pass(11) = isa(tubeTech(), 'trigtech') ; 

% Test making a chebfun3 from a scalar coefficient:
f = chebfun3(1, 'coeffs'); 
pass(12) = norm(f - 1) < tol ;

% Test passing an 'eps' value. Make sure it affects both rank and lengths.
ff = @(x,y,z) -x.*sin(sqrt(x+15+2*y+3*z));
f = chebfun3(ff);
fEps = chebfun3(ff, 'eps', 1e-8);
[m, n, p] = length(f);
[mEps, nEps, pEps] = length(fEps);
pass(13) = mEps < m && nEps < n && pEps < p ;

[r1, r2, r3] = rank(f);
[r1Eps, r2Eps, r3Eps] = rank(fEps);
pass(14) = r1Eps < r1 && r2Eps < r2 && r3Eps < r3 ;

% Construct from a string of constant type:
f = chebfun3('pi');
pass(15) = f(0,0,0) - pi < tol;

% Construct from a string of a univariate function
f = chebfun3('cos(alpha)');
pass(16) = f(0,0,0) - cos(0) < tol;

% Construct from a string of a bivariate function
f = chebfun3('x+y');
pass(17) = f(0.25,0.5,0) - 0.75 < tol;

% Construct from a string of a trivariate function
f = chebfun3('cos(x+y+z)');
pass(18) = f(0.25,0.5,-0.3) - cos(0.25+0.5-0.3) < tol;

end