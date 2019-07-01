function pass = test_sum3( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

% Example 1
f = ballfun(@(r,lam,th)1, 'spherical');
I = sum3(f);
exact = 4*pi/3;
pass(1) = abs(I-exact)<tol;

% Example 2
f = ballfun(@(x,y,z)x+1);
I = sum3(f);
exact = 4*pi/3;
pass(2) = abs(I-exact)<tol;

% Example 3
f = ballfun(@(x,y,z)x.^2);
I = sum3(f);
exact = 4*pi/15;
pass(3) = abs(I-exact)<tol;

if (nargout > 0)
    pass = all(pass(:));
end

end
