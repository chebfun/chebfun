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

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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

numCols = size(A,2);

if ( numel(A) > 1 )
    % Quasimatrix case:

    isSimple = true;
    for k = 1:numel(A)
        %isSimple = isSimple && all(cellfun(@(f) isa(f.onefun, 'chebtech'), A(k).funs));
        isSimple = isSimple && ~isdelta(A(k)) ...
            && ( all(cellfun(@(f) isa(f.onefun, 'chebtech'), A(k).funs))...
            ||   all(cellfun(@(f) isa(f.onefun, 'trigtech'), A(k).funs)) );
    end

    if ( isSimple )
        % Simple case = all FUNs are simple (i.e., use CHEBTECHs).
        % Convert to an array-valued CHEBFUN:
        A = quasi2cheb(A);
        [Q, R] = qr(A, 0);
        
    else
        % Legendre matrix:
        E = legpoly(0:numCols-1, domain(A), 'norm', 1);
        E = restrict(E, get(A, 'domain'));
        % Convert the Legendre-Vandermonde matrix to a quasimatrix:
        E = cheb2quasi(E);
        % Call abstract QR:
        [Q, R] = abstractQR(A, E, @innerProduct, @normest);
    end

elseif ( numel(A.funs) == 1 )
    % No breakpoints = easy case.
    
    % Call QR at the FUN level:
    [Q, R] = qr(A.funs{1});
    Q = chebfun({Q});
    
elseif ( size(A, 2) == 1 )
    % If f has only one column we simply scale it.
    R = sqrt(innerProduct(A, A));
    Q = A./R;

elseif ( all(cellfun(@(f) isa(f.onefun, 'chebtech'), A.funs)) )
    % Simple case = all FUNs are simple (i.e., use CHEBTECHs).
    %   (If all the FUN objects have .onefuns which are CHEBTECHs, we can use a
    %   much more efficient approach. Note, this completely violates OOP
    %   principles, but the performance gain is worth it.)

    % NOTE: We must use STRCMP here to get the expected behaviour.
    if ( strcmp(class(A.funs{1}.onefun), 'chebtech1') )  %#ok<STISA>
        chebType = 1;
    else
        chebType = 2;
    end
    [Q, R] = qrSimple(A, chebType);
    
else
    % Work in the general continuous setting:
    
    % Legendre-Vandermonde matrix:
    E = legpoly(0:numCols-1, A.domain, 'norm', 1);
    [A, E] = overlap(A, E);
    % Call abstract QR:
    [Q, R] = abstractQR(A, E, @innerProduct, @normest);
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, R] = qrSimple(A, chebType)
% This implementation is fast, but relies on the fact that everything is a
% CHEBTECH on a bounded domain at heart.
A = simplify(A);
    
% Get some useful values
numCols = numColumns(A);
tol = epslevel(A)*vscale(A);
dom = A.domain;
a = dom(1);
b = dom(end);
numFuns = length(dom) - 1;

% Get the sizes of the funs in the columns of A, keeping in mind that we
% will have to multiply with the columns of E and the norm of A's columns.
sizes = zeros(numFuns, 1);
for j = 1:numFuns
    sizes(j) = 2*max(length(A.funs{j}), numCols);
    A.funs{j}.onefun = prolong(A.funs{j}.onefun, sizes(j));
end

% Create the Chebyshev nodes and quadrature weights:
[pts, w] = chebpts(sizes, dom, chebType);

% Define the inner product as an anonymous function:
ip = @(f, g) w * (conj(f) .* g);

% Make the discrete analog of A:
A = get(A, 'values');
if ( iscell(A) )
    A = cell2mat(A);
end

% Generate a discrete E (Legendre-Chebyshev-Vandermonde matrix) directly:
xx = (pts - a)/(b - a) - (b - pts)/(b - a);  % Unscale the Chebyshev points.
E = ones(size(A));
E(:,2) = xx;
for k = 3:numCols % Recurrence relation:
    E(:,k) = ((2*k - 3)*xx.*E(:,k-1) - (k - 2)*E(:,k-2)) / (k - 1);
end
% Scaling:
for k = 1:numCols
    E(:,k) = E(:,k) * sqrt((2*k - 1) / (b - a));
end

[Q, R] = abstractQR(A, E, ip, @(v) norm(v, inf), tol);

% Construct a CHEBFUN from the discrete values:
pref = chebfunpref();
if ( chebType == 1 )
    pref.tech = @chebtech1;
else
    pref.tech = @chebtech2;
end
Q = mat2cell(Q, sizes, numCols);
Q = chebfun(Q, dom, pref);

end
