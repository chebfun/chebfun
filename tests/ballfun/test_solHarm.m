function pass = test_solHarm( pref ) 
% Test the equality for P^m_l, 0 <= l <= n and -m <= l <= m

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e4*pref.techPrefs.chebfuneps;

n = 10;
Max_difference = 0;

for l = 0:n
    for m = 0:l
        Y = extract_spherefun(ballfun.solHarm(l,m))*sqrt(2*l+1)*sqrt(1+(m>0));
        Z = spherefun.sphharm(l,m);
        Max_difference = max(norm(Y-Z),Max_difference);
    end
end
pass = Max_difference < tol;
end