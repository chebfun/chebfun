function pass = test_HelmholtzSolver( ) 
% Test the Helmholtz Solver on the sphere 

tol = 1e-10; 
for k = [1 2 3 4]
    
    if k == 1
        m = 60; n = 40;
    elseif k == 2
        m = 61; n = 41;
    elseif k == 3
        m = 62; n = 42;
    else
        m = 63; n = 43;
    end

    nxt = 1;  
    for L = 0:3
        for M=0:L 

        f = spherefun.sphharm(L, M); 

        K = 100.1; 
        u = spherefun.helmholtz((K^2-L*(L+1))*f, K, m, n);

        pass(nxt,k) = ( norm( u - f, 2 ) < 100*tol );

        nxt = nxt + 1; 
        end
    end

end

if ( nargout > 0 )
    pass = all(pass(:)); 
end

end