function g = squeeze(f)
%SQUEEZE  squeeze a chebfun2 to one variable, if possible.
% 
% G = squeeze(F) returns a chebfun2 if F depends on x and y. If F depends
% only on the x-variable a row chebfun is returned and if it depends on
% just the y-variable a column chebfun is returned. 

if ( isempty( f ) )  % check for an empty chebfun2.
   g = chebfun2;  
   return
end

if ( rank( f ) == 1 )   % must be of rank 1. 
    if ( norm( f.cols - mean(f.cols)) < 10*eps )
        g = mean(f.cols) * diag(1./f,pivotValues) * f.rows; 
        id = 1; 
    elseif ( norm( f.rows - mean(f.rows)) < 10*eps )
        g = f.cols * diag(1./f.pivotValues) * mean(f.rows);
        id = 3;
    else
        g = f; 
    end
    if ( isa(g,'double') )
        g = chebfun( g, [ rect(id), rect(id+1) ] );
    end
else
    g = f; 
end

if ( isa( g, 'chebfun' ) )
    g = simplify( g );
end

end