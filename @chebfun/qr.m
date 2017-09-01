function [Q, R] = qr(A, econ)
%QR   QR factorization of an array-valued CHEBFUN.
%   [Q, R] = QR(A) or QR(A, 0), where A is a column CHEBFUN with n columns,
%   produces a column CHEBFUN Q with n orthonormal columns and an n x n upper
%   triangular matrix R such that A = Q*R.
%
% See also SVD, MRDIVIDE, RANK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check inputs:
if ( (nargin == 2) && (econ ~= 0) )
    error('CHEBFUN:CHEBFUN:qr:twoargs',...
      ['Use qr(A) or qr(A, 0) for QR decomposition of an array-valued ' ...
      'CHEBFUN or quasimatrix']);
elseif ( A(1).isTransposed )
    error('CHEBFUN:CHEBFUN:qr:transpose',...
        'CHEBFUN QR works only for column CHEBFUN objects.')
elseif ( ~all(isfinite(A(1).domain)) )
    error('CHEBFUN:CHEBFUN:qr:infdomain', ...
        'CHEBFUN QR does not support unbounded domains.');
end

numCols = numColumns(A);
if ( numCols == 1 )
    % Trivial case: If A has only one column we simply scale it.
    R = sqrt(innerProduct(A, A));
    if ( R ~= 0 )
       Q = A./R;
    else
       Q = 1./sqrt(diff(A.domain)) + 0*A;
    end
    return
end
    
% Attempt to convert to an array-valued CHEBFUN:
[A, isArrayValued] = quasi2cheb(A);
    
if ( isArrayValued && (numel(A.funs) == 1) )
    % Array-valued CHEBFUN with a single FUN.
    
    % Call QR at the FUN level:
    [Q, R] = qr(A.funs{1});
    Q = chebfun({Q});

elseif ( isArrayValued )   
    % Array-valued CHEBFUN with multiple FUNS.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Developer note:
    %   Here we essentially use a panel-factored QR which allows us to do a QR
    %   factorization on each fun individually and then combine the result.
    %   Here's an example of this in a 2-FUN case:
    %    [A1] = [Q1*R^1] = [Q1 0][R^1] = [Q1 0][Q^1 ~][R] = [Q1*Q^1]R
    %    [A2] 1 [Q2*R^2]   [0 Q2][R^2] 2 [0 Q2][Q^2 ~][0] 3 [Q2*Q^2]
    %                                  ^
    %                               here [Q^:=Qhat, R] = qr(Rhat:=[R^1;R^2])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numFuns = numel(A.funs);

    % Step 1: Perform QR on each piece.
    Q = cell(numFuns, 1); Rhat = Q;
    for k = 1:numFuns
        [Q{k}, Rhat{k}] = qr(A.funs{k});
    end
    
    % Step 2: Compute [Qhat, R] = qr(Rhat),
    [Qhat, R] = qr(cell2mat(Rhat));
    R = R(1:numCols,:);       % Extract first block row.
    Qhat = Qhat(:,1:numCols); % Extract first block column.

    % Step 2b: Ensure the diagonal is non-negative. (A = QR = (Q*S)*(S*R))
    s = sign(diag(R)); s(~s) = 1;
    S = spdiags(s, 0, numCols, numCols);
    Qhat = Qhat*S;
    R = S*R;

    % Step 2c: Separate the segments of Qhat back into a cell.
    m = cellfun(@(v) size(v, 1), Rhat); % m(k) = length of A.FUN{k}.
    Qhat = mat2cell(Qhat, m, numCols);
    
    % Step 3: Fold Qhat back in to Q.
    Q = cellfun(@mtimes, Q, Qhat, 'UniformOutput', false);
    
    % Construct a new CHEBFUN from the computed FUNS:
    Q = chebfun(Q);
    
else
    % Quasimatrix case (tricky/slow):

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Developer note:
    %   Currently (5th Feb 2016) the only way this is reachable is in the
    %   case of a quasimatrix consisting of SINGFUN or DELTAFUN objects,
    %   neither of which return anything sensible when we attempt to compute
    %   a QR factorization. Try for example
    %    x = chebfun('x', [0 1]);
    %    qr([1 x sqrt(x)])
    %
    %   This case is not tested (which is OK)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Legendre-Vandermonde matrix:
    L = legpoly(0:numCols-1, domain(A), 'norm', 1);
    % Convert so that L is also quasimatrix:
    L = cheb2quasi(L);
    % Call abstract QR:
    [Q, R] = abstractQR(A, L, @innerProduct, @normest);

end

end
