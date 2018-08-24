function pass = test_cos( pref ) 
S = [20,21,22];

% Example 1
f = cos(ballfun(@(r,lam,th)r,S));
exact = ballfun(@(r,lam,th)cos(r),S);
pass(1) = norm( f - exact ) < tol ;

if (nargout > 0)
    pass = all(pass(:));
end
end
