function pass = test_sphharm( ) 
% Test spherical harmonic command.

tol = 1e4*chebfunpref().techPrefs.eps;

jj = 1; 
% Low order Y's:
for m = 1:10 
    for ll = m:m+3
        u = spherefun.sphharm( ll, m );
        v = laplacian( u ); 
        pass(jj) = ( norm( v + ll*(ll+1)*u, inf ) < tol );
        jj = jj + 1; 
    end
end 

end