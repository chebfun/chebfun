function pass = test_SeparableFormat( pref ) 
% Test code for computing separable format of PDOs. 
%
% Alex Townsend, September 2014. 

% Constant coefficient, Laplace: 
N = chebop2(@(x,y,u) diff(u, 2, 1) + diff(u, 2, 2), [-1 1 -1 1]); 
[U, S, V] = chebop2.separableFormat( N );
pass(1)  = checkCoeffs(N.coeffs, U, S, V); 

% Constant coefficient, helmholtz: 
N = chebop2(@(x,y,u) diff(u, 2, 1) + diff(u, 2, 2) + 100*u, [-1 1 -1 1]); 
[U, S, V] = chebop2.separableFormat( N );
pass(2)  = checkCoeffs(N.coeffs, U, S, V);

% Rank 3 PDO: 
N = chebop2(@(x,y,u) diff(u, 2, 1) + diff(u, 2, 2) + (10*x + y).*u, [-1 1 -1 1]); 
[U, S, V] = chebop2.separableFormat( N );
pass(3)  = checkCoeffs(N.coeffs, U, S, V);

% Highly oscillatory variable coefficient: 
N = chebop2(@(x,y,u) diff(u, 2, 1) + diff(u, 2, 2) + cos(100*x).*u, [-1 1 -1 1]); 
[U, S, V] = chebop2.separableFormat( N );
pass(4)  = checkCoeffs(N.coeffs, U, S, V);

% Helmholtz with gravity: 
N = chebop2(@(x,y,u) diff(u, 2, 1) + diff(u, 2, 2) + (10 + y).*u, [-1 1 -1 1]); 
[U, S, V] = chebop2.separableFormat( N );
pass(5)  = checkCoeffs(N.coeffs, U, S, V);

% On non-standard domain: 
N = chebop2(@(x,y,u) diff(u,2,1) + diff(u,2,2) + (10+y).*u, ...
    [-1 pi/3 -sqrt(3)/2 1.1]); 
[U, S, V] = chebop2.separableFormat( N );
pass(6)  = checkCoeffs(N.coeffs, U, S, V);

% Variable coefficients on higher derivatives: 
N = chebop2(@(x,y,u) (1+x.^2).*diff(u,2,1) + diff(u,2,2) + u, ...
    [-1 pi/3 -sqrt(3)/2 1.1]); 
[U, S, V] = chebop2.separableFormat( N );
pass(7)  = checkCoeffs(N.coeffs, U, S, V);

% Heat equation: 
N = chebop2(@(x,y,u) diff(u,1,1) - diff(u, 2, 2), [-1 1 0 10]); 
[U, S, V] = chebop2.separableFormat( N );
pass(8) = checkCoeffs(N.coeffs, U, S, V);

% Wave equation: 
N = chebop2(@(x,y,u) diff(u, 1, 2) - x.^2.*diff(u, 2, 2) + y.*u, [-1 1 0 10]); 
[U, S, V] = chebop2.separableFormat( N );
pass(9) = checkCoeffs(N.coeffs, U, S, V);

% complex scalar: 
N = chebop2(@(x,y,u) 1i*diff(u, 1, 1)); 
[U, S, V] = chebop2.separableFormat( N );
pass(10) = checkCoeffs(N.coeffs, U, S, V);

% Complex heat equation? 
N = chebop2(@(x,y,u) 1i*diff(u,1,1) + (1+1i)*diff(u,1,2)); 
[U, S, V] = chebop2.separableFormat( N );
pass(11) = checkCoeffs(N.coeffs, U, S, V);

% Time-independent Schrodinger (complex coefficients): 
hb = 0.0256; 
N = chebop2(@(x,y,u) 1i*hb*diff(u, 1, 1) + hb^2*diff(u, 2, 2) - x.^2.*u, ...
    [-1 1 0 10]); 
[U, S, V] = chebop2.separableFormat( N );
pass(12) = checkCoeffs(N.coeffs, U, S, V);

end

function bol = checkCoeffs( A, cellU, S, cellV )
% CheckCoeffs  debugging script. 
% Alex Townsend, September 2014. 
for jj = 1:size(cellU, 1)
    for kk = 1:size(cellV, 1)
        a = cellU{jj,1} * S(1,1) * cellV{kk,1}.';
        for r = 2:size(cellU,2)
            a = a + cellU{jj,r} * S(r,r) * cellV{kk,r}.';
        end
        store{jj,kk} = a; 
        err(jj,kk) = norm( a - A{jj,kk} );
    end
end
if ( norm( err ) < 1e-8 ) 
    bol = 1; 
else 
    bol = err; 
end
% debug plot: 
% sp = 1; 
% M = size(A, 1); 
% N = size(A, 2); 
% for jj = 1:M
%     for kk = 1:N
%         subplot(M,N,sp), plot(store{jj,kk} )
%         sp = sp + 1; 
%     end
% end
end

    

