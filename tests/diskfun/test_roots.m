function pass = test_roots( pref ) 
% Check that roots works for a diskfun.

if ( nargin < 1 ) 
    pref = chebfunpref; 
end
tol = 1e5*pref.cheb2Prefs.chebfun2eps;

% More extensive testing is needed.  Here we just test that errors aren't
% thrown and the form of the output is reasonable.

f = diskfun(@(x,y) x+y );
r = roots(f);
if ~isempty(r)
    pass(1) = 1;
else
    pass(1) = 0;
end

% No zero contour exists
f = diskfun(@(x,y) 2+x.^2);
r = roots(f);
if isempty(r)
    pass(2) = 1;
else
    pass(2) = 0;
end

% Zero contour at  x.^2+y.^2=.5^2
f = diskfun(@(t,r) r.^2-.5^2, 'polar');
r = roots(f);
% First component of roots should be cos(pi*lam) (or its negative).
r1true = chebfun(@(t) .5*cos(pi*t));
pass(3) = ( ( norm(r{1}(:,1)-r1true) < tol ) || ...
            ( norm(r{1}(:,1)+r1true) < tol ) );
% Second component of roots should be sin(pi*lam) (or its negative).
r2true = chebfun(@(t) .5*sin(pi*t));
pass(4) = ( ( norm(r{1}(:,2)-r2true) < tol ) || ...
            ( norm(r{1}(:,2)+r2true) < tol ) );

end
