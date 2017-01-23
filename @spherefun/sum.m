function f = sum( f, dim )
%SUM   Definite Integration of a SPHEREFUN.
%   G = sum(F,DIM) where DIM is 1 or 2 integrates only over THETA (latitude)
%   or LAMBDA (longitude) respectively,
%   and returns as its output a chebfun in the remaining variable.
%
%   G = sum(F) is the same as sum(F,1)
%
% See also SUM2. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    f = []; 
    return; 
end

% Default to y direction: 
if ( nargin == 1 )
    dim = 1;
end

% Get the low rank representation for f. 
[cols, D, rows] = cdr(f);
dom = f.domain; 

if ( dim == 1 )
    % Integration over y requires including the measure of the sphere.
    % Additionally, we should only integrate over half the domain of the
    % columns since the cols are doubled up.  The code below does these
    % integrals in a fast way.  See the comment in spherefun/sum2 for why
    % we do it this way.
    [a, ignore] = trigcoeffs(cols);
    k = (0:size(a, 1)-1).';
    intFactor = 2./(1-k(1:2:end).^2);
    intCols = sum(bsxfun(@times, a(1:2:end, :), intFactor));
    f = rows * (intCols * D).';
    if ( isa(f, 'chebfun') ) 
        f = simplify( f.', [], 'globaltol' ); 
    else
        % f = double 
        f = chebfun(f, dom(1:2)).'; 
    end
elseif ( dim == 2 )
    f = cols * ( D * sum( rows ).' );
    if  ( isa(f, 'chebfun') ) 
        f = simplify( f, [], 'globaltol' );
        % Restrict the domain to remove the doubling up in theta
        f = restrict( f, [0 pi]);
    else
        % f = double 
        f = chebfun( f, [0 pi] ); 
    end
else 
    error('CHEBFUN:SPHEREFUN:sum:unknown', ...
          'Undefined function ''sum'' for that dimension');
end

end
