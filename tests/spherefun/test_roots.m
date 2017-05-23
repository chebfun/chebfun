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
% Third component of roots should be zero.
pass(3) = ( norm(r{1}(:,3)) < tol );

% Zero contour along the y-axis
f = spherefun(@(x,y,z) x);
r = roots(f);
% First component of roots should be zero.
pass(4) = ( norm(r{1}(:,1)) < tol );

% Zero contour across the x-axis
f = spherefun(@(x,y,z) y);
r = roots(f);
% First component of roots should be zero.
pass(5) = ( norm(r{1}(:,2)) < tol );

% Check that the roots are actually zero curves of f.
f = spherefun(@(x,y,z) 2*sinh(5*x.*y.*z)); r = roots(f);
for j=1:length(r)
    t = sample(r{j});
    fvals = f(t(:,1),t(:,2),t(:,3));
    pass(5+j) = ( norm(fvals,inf) < 1e-3 );
end

end
