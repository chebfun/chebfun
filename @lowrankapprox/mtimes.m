function h = mtimes(f, g)
%*	CHEBFUN2 multiplication.
%   c*F or F*c multiplies a CHEBFUN2 F by a scalar c.
%
%   F*G computes the integral of F(s,y)G(x,s) over s, and this is the continuous
%   analogue of matrix-matrix multiplication.
%
% See also TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'chebfun2') )           % CHEBFUN2 * ???
    
    if ( isa(g, 'double') )         % CHEBFUN2 * DOUBLE
        if ( numel(g) == 1 )
            h = f;
            h.pivotValues = h.pivotValues ./ g;
        else
            error('CHEBFUN:CHEBFUN2:mtimes:size', 'Sizes are inconsistent.');
        end
        
    elseif ( isa(g, 'chebfun') )    % CHEBFUN2 * CHEBFUN
        cols = f.cols;
        rows = f.rows;
        fScl = diag( 1./f.pivotValues );
        X = innerProduct( rows, g );
        h = cols * fScl * X;
        
    elseif ( isa(g, 'chebfun2') )   % CHEBFUN2 * CHEBFUN2  compute
        
        % Get the columns and rows of f
        fCols = f.cols;
        fRows = f.rows;
        fScl = diag( 1./f.pivotValues );
        
        % Get the columns and rows of g
        gCols = g.cols;
        gRows = g.rows;
        gScl = diag( 1./g.pivotValues );
        
        % Compute integral in s.
        X = innerProduct( fRows, gCols );
        
        % Construct low rank form of the result:
        [U, S, V] = svd( fScl * X * gScl );
        h = f;
        h.cols = fCols * U;
        h.rows = gRows * V;
        h.pivotValues = 1 ./ diag( S );
        
    elseif ( isa(g, 'chebfun2v') )  % CHEBFUN * CHEBFUN2V
        nG = g.nComponents; 
        h = g; 
        gc = g.components; 
        for jj = 1:nG 
           h.components{jj} = times(f, gc{jj}); 
        end
        
    else
        error('CHEBFUN:CHEBFUN2:mtimes:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
        
    end
    
else
    
    h = mtimes(g.', f.').';
    
end

end
