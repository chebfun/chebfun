function pass = test_roots( pref ) 
% Check that roots works for a spherefun.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 10*pref.cheb2Prefs.chebfun2eps;

% More extensive testing is needed.  Here we just test that errors aren't
% thrown and the form of the output is reasonable.

% Zero contour at the equator
f = spherefun(@(x,y,z) z);
r = roots(f);
if ~isempty(r)
    pass(1) = 1;
else
    pass(1) = 0;
end

% No zero contour exists
f = spherefun(@(x,y,z) 2+z);
r = roots(f);
if isempty(r)
    pass(2) = 1;
else
    pass(2) = 0;
end

end
