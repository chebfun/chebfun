function pass = test_constructor2( pref ) 
% This tests the chebfun2 constructor. 

if ( nargin < 1 ) 
    pref = chebfunpref; 
end 
tol = 1e2 * pref.cheb2Prefs.chebfun2eps;


% Test building Chebfun2 objects from sample data: 
seedRNG(0);
r = rand(3);
pass(1) = norm( r - chebcoeffs2(chebfun2(r, 'coeffs')) ) < 10*tol; 
r = rand(4);
pass(2) = norm( r - chebcoeffs2(chebfun2(r, 'coeffs')) ) < 10*tol;
r = rand(4);
pass(3) = norm( r - chebpolyval2(chebfun2(r)) ) < 10*tol;

g = @(x,y) exp(-((x+pi).^2+y.^2)./max(1 - ((x+pi).^2 + y.^2),0));
f = chebfun2(g,[-pi,pi,-pi,pi]);
[xx,yy] = meshgrid(linspace(-pi,pi,1001));
err = g(xx,yy)-f(xx,yy);
pass(4) = ( norm(err(:),inf ) < 2e7*tol );

% Make a chebfun2 based on TRIGTECH by calling the 'trig' flag:
f1 = chebfun2(@(x,y) cos(pi*x).*sin(pi*y),'trig'); 
f2 = chebfun2(@(x,y) cos(pi*x).*sin(pi*y),[-1 1 -1 1],'trig');
pass(5) = ( norm( f1 - f2 ) < tol );

f1 = chebfun2(@(x,y) cos(pi*cos(pi*x) + pi*sin(pi*y)),'trig');
f2 = chebfun2(@(x,y) cos(pi*cos(pi*x) + pi*sin(pi*y)), [-1 1 -1 1], 'trig');
pass(6) = ( norm(f1 - f2) < 10*tol );
% Check underlying tech is a TRIGTECH: 
techRow = get(f1.cols.funs{1}, 'tech');
techCol = get(f1.rows.funs{1}, 'tech');
pass(7) = ( isa(techRow(), 'trigtech') ); 
pass(8) = ( isa(techCol(), 'trigtech') ); 

% Make sure the 'periodic' flag works as well:
f1 = chebfun2(@(x,y) cos(pi*cos(pi*x) + pi*sin(pi*y)),'periodic');
f2 = chebfun2(@(x,y) cos(pi*cos(pi*x) + pi*sin(pi*y)), [-1 1 -1 1], 'periodic');
pass(9) = ( norm(f1 - f2) < 10*tol );
% Check underlying tech is a TRIGTECH:
techRow = get(f1.cols.funs{1}, 'tech');
techCol = get(f1.rows.funs{1}, 'tech');
pass(10) = ( isa(techRow(), 'trigtech') );
pass(11) = ( isa(techCol(), 'trigtech') );

% Test making a chebfun2 from a scalar coefficient: 
f = chebfun2(1,'coeffs'); 
pass(12) = ( norm( f - 1 ) < tol );
f = chebfun2('x'); z = .5 + sqrt(3)/3*1i;
pass(13) = ( norm( f(real(z),imag(z)) - z ) < tol ) ; 

% Test making a chebfun2 with a fixed length:
m = 8;
n = 10;
f = chebfun2( @(x,y) cos(x.*y), [m n]);
[mF, nF] = length(f);
pass(14) = ( mF == m && nF == n );

% Test making a chebfun2 with a fixed rank
r = 5;
f = chebfun2( @(x,y) exp(cos(x.*y)), r );
pass(15) = ( rank(f) == r );

% Test making a chebfun2 with fixed rank and length
r = 2; m = 8; n = 10;
f = chebfun2( @(x,y) exp(cos(x.*y)), r, [m n] );
[mF, nF] = length(f);
pass(16) = ( rank(f) == r && mF == m && nF == n );

% Reverse the previous calling sequence.
f = chebfun2( @(x,y) exp(cos(x.*y)), [m n], r );
[mF, nF] = length(f);
pass(17) = ( rank(f) == r && mF == m && nF == n );

% Test making a chebfun2 with fixed rank and domain
dom = [-1.5 1.5 -0.5 0.75];
r = 3;
f = chebfun2( @(x,y) exp(cos(x.*y)), dom, r );
pass(18) = ( rank(f) == r && all( f.domain == dom ) );

% Reverse the previous calling sequence.
f = chebfun2( @(x,y) exp(cos(x.*y)), r, dom);
pass(19) = ( rank(f) == r && all( f.domain == dom ) );

% Test making a chebfun2 with fixed domain, fixed length, and rank
dom = [-1.5 1.5 -0.5 0.75];
r = 3; m = 20; n = 37;
f = chebfun2( @(x,y) exp(cos(x.*y)), dom, [m n], r );
[mF, nF] = length(f);
pass(20) = ( rank(f) == r && all( f.domain == dom ) && ...
             mF == m && nF == n );

% Test construction from a chebfun2 with fixed rank and length gives
% correct results.
g = chebfun2( @(x,y) exp(cos(x.*y)) );
r = 2; m = 8; n = 10;
f = chebfun2( g, r, [m n] );
[mF, nF] = length(f);
pass(21) = ( rank(f) == r && mF == m && nF == n );

% Test construction from a chebfun2 with fixed length and domain
% correct results.
m = 8; n = 10;
f = chebfun2( @(x,y) exp(cos(x.*y)), [m n], dom );
[mF, nF] = length(f);
pass(22) = ( all( f.domain == dom ) && mF == m && nF == n );

% Reverse
m = 8; n = 10;
f = chebfun2( @(x,y) exp(cos(x.*y)), dom, [m n] );
[mF, nF] = length(f);
pass(23) = ( all( f.domain == dom ) && mF == m && nF == n );

end