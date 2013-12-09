function h = outerProduct(f, g)
% Compute the outerProduct of two chebfuns. 

h = chebfun2;  % empty chebfun2 object 

cols = f; 
rows = g; 

% check the sizes are correct. 
if ( size(cols, 2) == size(rows, 1) )
    % check that things have been compressed.
    [Qleft, Rleft] = qr( cols );  
    [Qright, Rright] = qr( rows.' ); 
    [U, S, V] = svd( Rleft * Rright.' );
    nRank = find(diag(S) > max(f.epslevel, g.epslevel), 1, 'last'); 
    h.cols = Qleft * U(:, 1:nRank);
    h.rows = Qright * V(:, 1:nRank);
    s = diag(S); 
    h.pivotValues = 1 ./ s(1:nRank);
    h.domain = [g.domain f.domain];
else
    error('CHEBFUN2:OUTERPRODUCT:SIZES', 'Sizes not consistent for outerproduct');
end
end

%end