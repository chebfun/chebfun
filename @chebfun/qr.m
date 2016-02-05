function [Q, R] = qr(A, econ)
%QR   QR factorization of an array-valued CHEBFUN.
%   [Q, R] = QR(A) or QR(A, 0), where A is a column CHEBFUN with n columns,
%   produces a column CHEBFUN Q with n orthonormal columns and an n x n upper
%   triangular matrix R such that A = Q*R.
%
%   The algorithm used is described in L.N. Trefethen, "Householder
%   triangularization of a quasimatrix", IMA J. Numer. Anal. (30), 887-897
%   (2010).
%
% See also SVD, MRDIVIDE, RANK.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer note:
%  If A contains only a single FUN, then FUN/QR is used directly. If A has
%  multiple pieces but each of these are simple CHEBTECH objects, then
%  QRSIMPLE() is called. This violates OOP principles, but is _much_ more
%  efficient.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check inputs
if ( (nargin == 2) && (econ ~= 0) )
    error('CHEBFUN:CHEBFUN:qr:twoargs',...
      'Use qr(A) or qr(A, 0) for QR decomposition of an array-valued CHEBFUN.');
end
if ( A(1).isTransposed )
    error('CHEBFUN:CHEBFUN:qr:transpose',...
        'CHEBFUN QR works only for column CHEBFUN objects.')
end
if ( ~all(isfinite(A(1).domain)) )
    error('CHEBFUN:CHEBFUN:qr:infdomain', ...
        'CHEBFUN QR does not support unbounded domains.');
end

numCols = numColumns(A);

if ( numCols == 1 )
    % If f has only one column we simply scale it.
    R = sqrt(innerProduct(A, A));
    Q = A./R;
    
elseif ( numel(A) > 1 )
    % Quasimatrix case:

    % Try to convert to an array-valued CHEBFUN:
    [A, isSimple] = quasi2cheb(A);
    
    if ( isSimple ) % Successfully converted to array-valued CHEBFUN.
        
        [Q, R] = qr(A, 0);
        
    else            % We can't convert to array valued. Tricky case.
        
        % TODO: This case is not tested!
        %  DEVELOPER NOTE:
        %   Currently (5th Feb 2016) the only way this is reachable is in the
        %   case of a quasimatrix consisting of SINGFUN or DELTAFUN objects,
        %   neither of which return anything sensible when we attempt to compute
        %   a QR factorization.

        % Legendre-Vandermonde matrix:
        L = legpoly(0:numCols-1, domain(A), 'norm', 1);
        % Convert so that L is also quasimatrix:
        L = cheb2quasi(L);
        % Call abstract QR:
        [Q, R] = abstractQR(A, L, @innerProduct, @normest);
        
    end
    
elseif ( numel(A.funs) == 1 )
    % No breakpoints = easy case.
    
    % Call QR at the FUN level:
    [Q, R] = qr(A.funs{1});
    Q = chebfun({Q});

else

    numFuns = numel(A.funs);

    % Perform QR on each piece:
    Q = A;
    R = cell(numFuns, 1);
    for k = 1:numFuns
        [Q.funs{k}, R{k}] = qr(Q.funs{k});
    end

    % Compute [Qhat, Rhat] = qr(Q):
    m = cellfun(@(v) size(v, 1), R);
    [Qhat, Rhat] = qr(cell2mat(R));
    Rhat = Rhat(1:numCols,:); % Extract require entries.
    Qhat = mat2cell(Qhat(:,1:numCols), m, numCols);

    % Ensure diagonal has positive sign ...
    s = sign(diag(Rhat)); s(~s) = 1; 
    S = spdiags(s, 0, numCols, numCols); 
    R = S*Rhat;
    % ... and fold Qhat back in to Q:
    for k = 1:numFuns
        Q.funs{k} = Q.funs{k}*(Qhat{k,1}*S);
    end

end

end
