function pass = test_roots( pref ) 
% Check that roots works for a spherefun.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e5*pref.cheb2Prefs.chebfun2eps;

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

% Zero contour at the equator
f = spherefun(@(x,y,z) z);
r = roots(f);
% First component of roots should be cos(pi*lam) (or it's negative).
r1true = chebfun(@(lam) cos(pi*lam));
pass(3) = ( ( norm(r{1}(:,1)-r1true) < tol ) || ...
            ( norm(r{1}(:,1)+r1true) < tol ) );
% Second component of roots should be sin(pi*lam) (or it's negative).
r2true = chebfun(@(lam) sin(pi*lam));
pass(4) = ( ( norm(r{1}(:,2)-r2true) < tol ) || ...
            ( norm(r{1}(:,2)+r2true) < tol ) );
% Third component of roots should be zero.
pass(5) = ( norm(r{1}(:,3)) < tol );

end
