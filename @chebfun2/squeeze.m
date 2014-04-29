function g = squeeze(f)
%SQUEEZE   Squeeze a CHEBFUN2 to one variable, if possible.
%   G = squeeze(F) returns a CHEBFUN2 if F depends on x and y. If F depends only
%   on the x-variable a row CHEBFUN is returned and if it depends on just the
%   y-variable a column CHEBFUN is returned.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) )  % Check for an empty CHEBFUN2.
    g = chebfun();    % Return an empty CHEBFUN because we are squeezing.
    return
end

% Get the low rank representation for f.
[cols, D, rows] = cdr(f);
dom = f.domain;

if ( rank( f ) == 1 )
    % If f is of rank 1, then it may be a function of just one variable:
    
    % Simplify rows and cols:
    cols = simplify( cols );
    rows = simplify( rows );
    if ( length( cols ) == 1 )
        % If cols are constant then function of x.
        g = mean( cols ) * D * rows.';
        newdomain = dom(1:2);
    elseif ( length( rows ) == 1 )
        g = cols * D * mean(rows);
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
