function h = plus(f, g) 

% Chebfun2 plus using compression. 

if ( ~isa(f, 'chebfun2') )  % ??? + CHEBFUN2
    
    h = plus(g, f); 
    return
    
elseif ( isempty(g) )       % CHEBFUN2 + []

    f = [];
    
elseif ( isnumeric( g ) )   % CHEBFUN2 + DOUBLE
      g = chebfun2( g ); 
      
      h = plus( f, g ); 
      
      return
    
elseif ( ~isa(g, 'chebfun2') ) % CHEBFUN2 + ???

    error('CHEBFUN2:plus:unknown', ...
          ['Undefined function ''plus'' for input arguments of type %s ' ...
           'and %s.'], class(f), class(g));
       
elseif ( isempty(f) )         % empty CHEBFUN2 + CHEBFUN2

    % Nothing to do. (Return empty CHEBFUN2 as output).

else                          % CHEBFUN2 + CHEBFUN2
    h = compression_plus(f, g); 
end

end


function h = compression_plus(f, g)
% Add chebfun2 objects together by compression algorithm. 

% If A = XY^T and B = WZ^T, then A + B = [X W]*[Y Z]^T, 
% [Qleft, Rleft] = qr([X W])
% [Qright, Rright] = qr([Y Z])
% A = Qleft * (Rleft * Rright') * Qright'
% [U, S, V] = svd( Rleft * Rright' )
% A = (Qleft * U) * S * (V' * Qright')     -> new low rank representation
    
    h = f; 
    
    fScl = diag( 1./f.pivotValues );
    gScl = diag( 1./g.pivotValues ); 
    cols = [f.cols g.cols];
    rows = [f.rows g.rows];
    
    [Qleft, Rleft] = qr( cols ); 
    
    [Qright, Rright] = qr( rows );
    
    Z = zeros(length(fScl),length(gScl));  
    [U, S, V] = svd( Rleft * [fScl Z; Z' gScl] * Rright.' ); 
    
    h.cols = Qleft  * U; 
    h.rows = Qright * V;
    h.pivotValues = 1./diag( S ); 
end