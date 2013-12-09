function varargout = chebpolyval2( f, m, n )
% Evaluate a chebfun2 on a m-by-n Chebyshev-tensor grid

cols = f.cols; 
rows = f.rows; 
piv = f.pivotValues; 

C = resize( chebpoly( cols )', m);
R = resize( chebpoly( rows )', n);

C = chebtech2.coeffs2vals( C );
R = chebtech2.coeffs2vals( R )';

if nargout <= 1
    varargout = {C * diag(1./piv) * R}; 
else
    varargout = {C diag(1./piv) R}; 
end
    

end

function X = resize( X, N )
% Resize the matrix to have length N.
    
[mX, nX] = size( X ); 

if ( mX > N ) 
    % truncate 
    X = X(end-N+1:end, :); 
else
    % pad
    X = [zeros(N-mX, nX); X]; 
end

end