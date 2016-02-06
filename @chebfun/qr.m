function [Q, R] = qr(A, econ)
%QR   QR factorization of an array-valued CHEBFUN.
%   [Q, R] = QR(A) or QR(A, 0), where A is a column CHEBFUN with n columns,
%   produces a column CHEBFUN Q with n orthonormal columns and an n x n upper
%   triangular matrix R such that A = Q*R.
%
% See also SVD, MRDIVIDE, RANK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
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
    Q = A./R;
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
    % Array-valued CHEBFUN with multiple funs.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Developer note:
    %   Here we essentially use a panel-factored QR which allows us to do a QR
    %   factorization on each fun individually and then combine the result.
    %   Here's an example of this in a 2-FUN case:
    %    [A1] = [Q1*R1] = [Q1 0][R1] = [Q1 0][Q^1 ~][R^] = [Q1*Q^1]R^
    %    [A2]   [Q2*R2]   [0 Q2][R2]   [0 Q2][Q^2 ~][0 ]   [Q2*Q^2]
    %                                ^
    %                               here [Q^:=Qhat, R^:=Rhat] = qr(R:=[R1;R2])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numFuns = numel(A.funs);

    % Perform QR on each piece:
    R = cell(numFuns, 1);
    Q = R;
    for k = 1:numFuns
        [Q{k}, R{k}] = qr(A.funs{k});
    end

    [Qhat, Rhat] = qr(cell2mat(R));  % Compute [Qhat, Rhat] = qr(R). 
    Rhat = Rhat(1:numCols,:);        % Extract require entries.
    m = cellfun(@(v) size(v, 1), R); % m(k) = length of A.FUN{k}.
    Qhat = mat2cell(Qhat(:,1:numCols), m, numCols);

    % Ensure diagonal has positive sign. (A = QR -> (Q*S)*(S*R))
    s = sign(diag(Rhat)); s(~s) = 1; 
    S = spdiags(s, 0, numCols, numCols); 
    R = S*Rhat; % (Qhat --> Qhat*S is performed below)
    
    % Fold Qhat back in to Q:
    for k = 1:numFuns
        Q{k} = Q{k}*(Qhat{k,1}*S);
    end
    
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
