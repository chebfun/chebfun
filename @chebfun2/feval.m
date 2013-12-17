function out = feval(f, x, y)
% FEVAL  evaluation of a chebfun2 at a point. 

if ( strcmpi(x, ':') && strcmpi(y, ':') )
    out = f;
elseif ( strcmpi(x, ':') && isnumeric( y ) )
    Cols = f.cols;
    Rows = f.rows;
    pivotValues = f.pivotValues;
    y = y(:);
    
    out = feval( Cols, y ) * diag( 1./pivotValues ) * Rows';
    out = simplify( out ); 
elseif ( isnumeric( x ) && strcmpi(y, ':') )
    Cols = f.cols;
    Rows = f.rows;
    pivotValues = f.pivotValues;
    x = x( : ); 
    
    out = Cols * diag( 1./pivotValues ) * feval( Rows, x )';
    out = simplify( out ); 
elseif ( isnumeric( x ) && isnumeric( y ) )

    sx = size(x);
    sy = size(y);
    
    if ( min(sx) > 1 && all(sx == sy) )
        if ( rank(x) == 1 && rank(y) == 1 )
            x = x(1,:);
            y = y(:,1);
        end
    end
    
    zCol = feval(f.cols, y(:));
    zRow = feval(f.rows, x(:));
    
    out = zCol*diag(1./f.pivotValues)*zRow.';
else
    error('CHEBFUN2:FEVAL:INPUTS','Unrecognized arguments for evaluation');
end

end