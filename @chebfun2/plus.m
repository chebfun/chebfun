function h = plus(f, g)
%+	  Plus for CHEBFUN2 objects.
%
% F + G adds F and G. F and G can be scalars or CHEBFUN2 objects.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( ~isa(f, 'chebfun2') )      % ??? + CHEBFUN2
    
    h = plus(g, f);
    
elseif ( isempty( g ) || isempty( f ))  % CHEBFUN2 + []
    
    % Return empty CHEBFUN2.
    h = chebfun2();
    
elseif ( isa( g, 'double' ) )   % CHEBFUN2 + DOUBLE
    
    % Convert g to a CHEBFUN2
    g = chebfun2( g, f.domain );
    h = plus( f, g );
    
elseif ( ~isa(g, 'chebfun2') )  % CHEBFUN2 + ???
    
    error( 'CHEBFUN2:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class(f), class(g) );
    
else                            % CHEBFUN2 + CHEBFUN2
    
    % Domain Check:
    if ( ~domainCheck(f, g) )
        error('CHEBFUN2:PLUS:DOMAIN', 'Inconsistent domains.');
    end
    
    % Check for zero CHEBFUN2 objects:
    if ( iszero( f ) )
        h = g;
    elseif ( iszero( g ) )
        h = f;
    else
        % Add together two nonzero CHEBFUN2 objects
        h = compression_plus(f, g);
        %h = chebfun2(@(x, y) feval(f, x, y) + feval(g, x, y), f.domain); 
    end 
    
    
    % Add the derivatives. Need to make sure they have the correct dimensions.
    fderiv = f.deriv;
    gderiv = g.deriv;
    
    % Obtain the dimensions of the derivatives of the inputs.
    [mf, nf] = size(fderiv);
    [mg, ng] = size(gderiv);
    
    % The final derivative will have the dimensions corresponding to the maximum
    % dimensions of the input derivatives. We create an all zero matrix of that
    % dimensions, then add the old derivative information to the bottom right of
    % that matrix
    newDeriv = zeros(max(mf,mg),max(nf,ng));
    fderivNew = newDeriv;
    gderivNew = newDeriv;
    
    % Replace entries of the matrices of the correct size with information from
    % the input variables.
    fderivNew(end-mf+1:end,end-nf+1:end) = fderiv;
    gderivNew(end-mg+1:end,end-ng+1:end) = gderiv;
    
    % Add the two derivative matrices together, and assign to the deriv field of
    % the output variable.
    h.deriv = fderivNew + gderivNew;
    
    
    
end

end

function h = compression_plus(f, g)
% Add CHEBFUN2 objects together by a compression algorithm.

% The algorithm is as follows:
% If A = XY^T and B = WZ^T, then A + B = [X W]*[Y Z]^T,
% [Qleft, Rleft] = qr([X W])
% [Qright, Rright] = qr([Y Z])
% A + B = Qleft * (Rleft * Rright') * Qright'
% [U, S, V] = svd( Rleft * Rright' )
% A + B = (Qleft * U) * S * (V' * Qright')     -> new low rank representation

% Hack: Ensure g has the smaller pivot values.
if ( norm(f.pivotValues, -inf) < norm(g.pivotValues, -inf) )
    % TODO: Understand why this works!
    h = compression_plus(g, f);
    return
end

fScl = diag( 1./f.pivotValues );
gScl = diag( 1./g.pivotValues );
cols = [f.cols, g.cols];
rows = [f.rows, g.rows];

[Qcols, Rcols] = qr( cols );
[Qrows, Rrows] = qr( rows );

Z = zeros( length(fScl), length(gScl) );
D = [fScl Z ; Z.' gScl];
[U, S, V] = svd( Rcols * D * Rrows.' );

h = f;
h.cols = Qcols * U;
h.rows = Qrows * V;
h.pivotValues = 1./diag( S );

end

