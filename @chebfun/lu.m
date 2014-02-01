% TODO: 

f = chebfun2(@(x,y) cos(x.*y) + x); 
C = f.cols; 
[Q, R] = qr( C ); 

L = C; 
U = zeros( eye( size(L, 2) ) ); 
% Do GE with partial pivoting: 
for j = 1 : size(Q, 2)
    % find maximum absolute entry: 
    [vals, pos] = minandmax( L(:, j) );
    [mx, idx] = max( abs( vals ) ); 
    pos = pos( idx );
    
    % Do one step of GE: 
    L = L - L(:, j) * L(pos, :) / mx; 
    
    % Store upper-triangular part: 
    U = U + L(pos, :); 
    
    % store pivot locations: 
    pivots( j ) = pos; 
end


roots(L)

