function f = simplify( f, tol, rank_flag )
% Simplify a SEPARABLEAPPROX
%
% F = SIMPLIFY( F ) compressed the representation of F to one that is
% numerically the same, but requires fewer parameters to store. This
% simplifies the bivariate polynomial degree of F, but not its rank.
%
% F = SIMPLIFY(F, TOL) does the same as SIMPLIFY( F ) but uses the scalar 
% TOL instead of the default simplification tolerance as the relative 
% threshold level for compression.
%
% F = SIMPLIFY(F, 'rank') compressed the rank of the representation for F 
% to one that is numerically the same. 
%
% F = SIMPLIFY(F, TOL, 'rank') does the same as SIMPLIFY(F, 'rank') but 
% uses the scalar TOL instead of the default simplification tolerance as 
% the relative threshold level for compression.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

simplifyRank = 0;   % do not simplify the rank by default.
if ( nargin < 2 )
    tol = [];
end

if ( nargin < 3 ) 
    if ( strcmpi(tol, 'rank') )
        simplifyRank = 1; 
        pref = chebfunpref; 
        tol = pref.cheb2Prefs.chebfun2eps; 
%     else
%        tol = []; 
    end
end

if ( nargin == 3 && strcmpi(rank_flag, 'rank')  )
    simplifyRank = 1;
end

if ( ~simplifyRank )
    % Note that we do not simplify the rank here because that would require
    % calling the SVD, which is expensive.
    
    % Simplify the column and row slices.
    f.cols = simplify( f.cols, tol, 'globaltol' );
    f.rows = simplify( f.rows, tol, 'globaltol' );
    
    % Ensure the left and right limits match the endpoints:
    f.cols = resetPointValues(f.cols);
    f.rows = resetPointValues(f.rows);
else
    % To simplify the rank we must call the SVD: 
    [U, S, V] = svd( f ); 
    s = diag(S);
    idx = find( s > tol, 1, 'last'); 
    f.cols = U(:, 1:idx); 
    % [TODO]: Pivot locations have very little meaning after this 
    % compression step. We do not update them. Should we? 
    f.pivotValues = 1./s(1:idx);
    f.rows = V(:, 1:idx); 
end

end
