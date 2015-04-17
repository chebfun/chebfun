function f = compress( f )
% COMPRESS     Compress the rank of a spherefun
%
%  F = COMPRESS( F ), use the SVD to compress the rank of a spherefun.



% DISCLAIMER: This is a first effort and does need checking 

tol = 50*eps;
C = f.Cols.values;
R = f.Rows.values;
D = f.BlockDiag;
for j = 1:size(D,1)/2
    
    ii = 2*j-1:2*j;
    % Grady's memo tells us this preserves the BMC structure:
    [UC, SC, VC] = svd( C(:,ii) );
    [UR, SR, VR] = svd( R(:,ii) );
    
    [U, S, V] = svd( SC * VC' * D(ii,ii) * VR * SR ,0 );
    
    % Push back into CDR structure:
    D(ii,ii) = S;
    C(:,ii) = UC * U;
    R(:,ii) = UR * V;
end

% COMPRESS: Remove zeros on diag of B
idx = find( diag( D ) > tol );
D = D( idx, idx );
C = C(:, idx);
R = R(:, idx);

% Now make a new spherefun:
f.Cols = trigtech( C );
f.Rows = trigtech( R );
f.BlockDiag = D;


end