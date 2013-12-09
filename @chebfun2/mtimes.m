function h = mtimes(f, g)
% MTIMES for Chebfun2 


if ( isa(f, 'chebfun2') )     % CHEBFUN2 * ???
    if ( isa(g, 'double') )   % CHEBFUN2 * DOUBLE
        if ( numel(g) == 1 )
            h = f; 
            h.pivotValues = g * h.pivotValues; 
        else
           error('CHEBFUN2:MTIMES','Sizes are inconsistent.'); 
        end
    elseif ( isa(g, 'chebfun') )  % CHEBFUN2 * CHEBFUN 
        cols = f.cols; 
        rows = f.rows; 
        fScl = diag( 1./f.pivotValues ); 
        X = innerProduct( rows, g );
        h = cols * fScl * X;
    elseif ( isa(g, 'chebfun2') )  % CHEBFUN2 * CHEBFUN2 
        h = f; 
        
        fCols = f.cols; 
        fRows = f.rows; 
        fScl = diag( 1./f.pivotValues ); 
        
        gCols = g.cols; 
        gRows = g.rows; 
        gScl = diag( 1./g.pivotValues );
        
        X = innerProduct( fRows, gCols ); 
        
        [U, S, V] = svd( fScl * X * gScl );
        
        h.cols = f.cols * U;
        h.rows = g.rows * V;
        h.pivotValues = 1 ./ diag( S ); 
        
    else
        
    error('CHEBFUN2:MTIMES:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type %s ' ...
           'and %s.'], class(f), class(g));
    end
else   
    h = transpose( mtimes(transpose( g ), transpose( f ) ) );
end

end
