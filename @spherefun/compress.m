function f = compress( f )
% COMPRESS     Compress the rank of a spherefun
%
%  F = COMPRESS( F ), use the SVD to compress the rank of a spherefun.

tol = 50*eps;

D = f.BlockDiag;

% Check to see if f has already been compressed.  If D is strictly
% diagonal then it has been compressed.
[ii,jj] = find(D);
if norm(ii-jj,inf) == 0
    return
end

C = f.cols.values;
R = f.rows.values;

for j = 1 : size(D, 1)/2
    
    % Look at 2x2 blocks at one time. 2x2 blocks separate. 
    ii = 2*j-1:2*j;
    
    % Follows equation (4) of Grady's memo.  Basic idea is to compute
    % eigenvalue decomposition of each 2-by-2 matrix on the block diagonal
    % M = U*S*U' 
    %   = 1/sqrt(2)[1 1;1 -1]*[M(1,1)+M(1,2) 0;0 M(1,1)-M(1,2)]*1/sqrt(2)*[1 1;1 -1]
    % Then Compute C*U and U'*R.  This will preserve the BMC structure
    
    % TODO: Perhaps we shoudl switch to the SVD of M so all pivots are
    % positive?
    
    C(:,ii) = [C(:,ii(1))+C(:,ii(2)) C(:,ii(1))-C(:,ii(2))]/sqrt(2);
    R(:,ii) = [R(:,ii(1))+R(:,ii(2)) R(:,ii(1))-R(:,ii(2))]/sqrt(2);    
    D(ii,ii) = spdiags([D(ii(1),ii(1))+D(ii(1),ii(2));D(ii(1),ii(1))-D(ii(1),ii(2))],0,2,2);
        
end

% COMPRESS: Remove zeros on diag of B
D = diag( D );
idx = find( abs( D ) > tol );
D = D( idx );
C = C(:, idx);
R = R(:, idx);

% order the columns and rows:
[ignored, perm] = sort( abs( D ), 1, 'descend'); 
D = D(perm); 
C = C(:, perm);
R = R(:, perm); 

% Now make a new spherefun:
f.cols = trigtech( C );
f.rows = trigtech( R );
f.blockDiag = spdiags(D,0,size(D,1),size(D,1));

end