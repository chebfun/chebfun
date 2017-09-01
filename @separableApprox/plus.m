function h = plus(f, g)
%+   Plus for SEPARABLEAPPROX objects.
%
% F + G adds F and G. F and G can be scalars or SEPARABLEAPPROX objects.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(f, 'separableApprox') ) % ??? + SEPARABLEAPPROX
    
    h = plus(g, f);
    
elseif ( isempty(g) ) % SEPARABLEAPPROX + []
    
    h = f; 
    
elseif ( isempty(f) ) % [] + SEPARABLEAPPROX
    
    h = g; 
    
elseif ( isa( g, 'double' ) )           % SEPARABLEAPPROX + DOUBLE
    
    g = compose( 0*f,@plus, g);   % promote double to object class.  
    h = plus(f, g); 
    
elseif ( ~isa(g, 'separableApprox') )          % SEPARABLEAPPROX + ???
    
    error( 'CHEBFUN:SEPARABLEAPPROX:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class(f), class(g));
    
else                                     % SEPARABLEAPPROX + SEPARABLEAPPROX
    
    % Domain Check:
    if ( ~domainCheck(f, g) )
        error('CHEBFUN:SEPARABLEAPPROX:plus:domain', 'Inconsistent domains.');
    end
    
    % Type check:
    if ( ~strcmp( class(f),class(g) ) )
        error( 'CHEBFUN:SEPARABLEAPPROX:plus:unknown', ...
            ['Undefined function ''plus'' for input arguments of type %s ' ...
            'and %s. Try converting to the same type.'], class(f), class(g));
    end 
    
    % Check for zero SEPARABLEAPPROX objects:
    if ( iszero(f) )
        h = g;
    elseif ( iszero(g) )
        h = f;
    else
        % Add together two nonzero SEPARABLEAPPROX objects:
        h = compression_plus(f, g);

    end 
    
end

end

function h = compression_plus(f, g)
% Add SEPARABLEAPPROX objects together by a compression algorithm.

% The algorithm is as follows:
% If A = XY^T and B = WZ^T, then A + B = [X W]*[Y Z]^T,
% [Qleft, Rleft] = qr([X W])
% [Qright, Rright] = qr([Y Z])
% A + B = Qleft * (Rleft * Rright') * Qright'
% [U, S, V] = svd( Rleft * Rright' )
% A + B = (Qleft * U) * S * (V' * Qright')     -> new low rank representation

% Hack: Ensure g has the smaller pivot values.
if ( norm(f.pivotValues, -inf) < norm(g.pivotValues, -inf) )
    % [TODO]: Understand why this works!
    h = compression_plus(g, f);
    return
end

fScl = diag(1./f.pivotValues);
gScl = diag(1./g.pivotValues);
cols = [f.cols, g.cols];
rows = [f.rows, g.rows];

[Qcols, Rcols] = qr(cols);
[Qrows, Rrows] = qr(rows);

Z = zeros(length(fScl), length(gScl));
D = [ fScl, Z ; Z.', gScl ];
[U, S, V] = svd(Rcols * D * Rrows.');
% Take diagonal from SIGMA:
s = diag(S);

% Compress the format if possible.
% [TODO]: What should EPS be in the tolerance check below?
vf = vscale(f); 
vg = vscale(g);

vscl = 2*max(vf, vg); 
% Remove singular values that fall below eps*vscale: 
idx = find( s > 10*eps * vscl, 1, 'last');

if ( isempty(idx) )
    % Return 0 separableApprox
    h = 0*f;
else
    U = U(:,1:idx);
    V = V(:,1:idx);
    s = s(1:idx);
    h = f;
    h.cols = Qcols * U;
    h.rows = Qrows * conj(V);
    % [TODO]: PivotValues have very little meaning after this compression step.
    % For now we assign the singular values as the pivot values. 
    h.pivotValues = 1./s;
end

end
