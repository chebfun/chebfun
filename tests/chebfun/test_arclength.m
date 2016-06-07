% Test file for @chebfun/arclength.m.

function pass = test_arclength(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

% A simple example:
f = chebfun(@(x)sin(50*x));

l = arcLength(f);
% The exact result is obtained using method 'integral' of Matlab:
lExact = 63.549358387599270; 
err = abs(l - lExact);
tol = 1e6*lExact*vscale(f)*eps;
pass(1) = ( err < tol );

dom = [-0.5 0.7];
lPart = arcLength(f{dom(1), dom(2)});
% The exact result is obtained using method 'integral' of Matlab:
lPartExact = 38.339961332457015;
err = abs(lPart - lPartExact);
tol = 1e7*lPartExact*vscale(f)*eps;
pass(2) = ( err < tol );


% A piecewise smooth CHEBFUN:
f = chebfun({@(x)sin(x) @(x)cos(2*x)}, [-1 1 2]);
l = arcLength(f);
% The exact result is obtained using method 'integral' of Matlab:
lExact = 4.044985856867475; 
err = abs(l - lExact);
pass(3) = ( err < lExact*vscale(f)*eps );

% A complex-valued CHEBFUN - A unit circle in the complex plane:
f = chebfun(@(x) exp(2*pi*1i*x), [0 1]);
l = arcLength(f);
lExact = 2*pi;
err = abs(l - lExact);
pass(4) = ( err < 10*lExact*vscale(f)*eps );

% Array-valued CHEBFUN:
f = chebfun(@(x)[sin(x) cos(2*x)]);
l = arcLength(f, [-0.6 0.9]);
% The exact result is obtained using method 'integral' of Matlab:
lExact = [ 2.019977962648605 2.482606372493216];
err = abs(l - lExact);
pass(5) = ( norm(err, inf) < 1e3*max(lExact*vscale(f)*eps) );

% Test on finite SINGFUN:
f = chebfun(@(x) sin(2*x).*((x+1).^0.5), 'exps', [0.5 0], 'splitting', 'on');
l = arcLength(f);
% The exact result is obtained using method 'integral' of Matlab:
lExact = 3.452674964506957;
err = abs(l - lExact);
pass(6) = ( err < 1e3*lExact*vscale(f)*eps );


end
