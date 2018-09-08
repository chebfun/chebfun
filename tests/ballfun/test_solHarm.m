function pass = test_solharm( pref ) 
% Test the equality for P^m_l, 0 <= l <= n and -m <= l <= m

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end
tol = 1e2*pref.techPrefs.chebfuneps;

n = 10;
Max_difference = 0;
Max_norm = 0;

for l = 0:n
    for m = 0:l
        Yml = ballfun.solharm(l,m);
        Y = Yml(1,:,:)*sqrt(1+(m>0))/sqrt(2*l+3);
        Z = spherefun.sphharm(l,m);
        % Check that the Legendre polynomials is the same as in Spherefun
        Max_difference = max(norm(Y-Z),Max_difference);
        % Check that the 2-norm is 1
        Max_norm = max(abs(norm(Yml)-1),Max_norm);
    end
end
pass(1) = Max_difference < tol;
pass(2) = Max_norm < tol;

if (nargout > 0)
    pass = all(pass(:));
end
end