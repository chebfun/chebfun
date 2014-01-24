function g = squeeze(f)
%SQUEEZE  squeeze a chebfun2 to one variable, if possible.
% 
% G = squeeze(F) returns a chebfun2 if F depends on x and y. If F depends
% only on the x-variable a row chebfun is returned and if it depends on
% just the y-variable a column chebfun is returned. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) )  % check for an empty chebfun2.
   g = chebfun;      % Return an empty chebfun because we are squeezing.
   return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 
piv = f.pivotValues; 
d = 1./piv; 
d(d==inf) = 0;  % set infinite values to zero.
dom = f.domain; 

% If f is of rank 1, then it may be a function of just one variable:
if ( rank( f ) == 1 )                % If of rank 1.  
    cols = simplify( cols ); 
    rows = simplify( rows ); 
    if ( length( cols ) == 1 )       % If cols are constant then function of x.
        g = mean( cols ) * diag( d ) * rows'; 
        newdomain = dom(1:2); 
    elseif ( length( rows ) == 1 )
        g = cols * diag( d ) * mean(rows);
        newdomain = dom(3:4); 
    else
        g = f; 
    end
    if ( isa(g,'double') )
        g = chebfun( g, newdomain );
    end
else
    g = f; 
end

if ( isa( g, 'chebfun' ) )
    g = simplify( g );
end

end