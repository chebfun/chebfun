function h = mtimes(f, g)
%*	   Pointwise multiplication for SEPARABLEAPPROX objects.
%
%   c*F or F*c multiplies a SEPARABLEAPPROX F by a scalar c.
%
%   F*G computes the integral of F(s,y)G(x,s) over s, and this is the continuous
%   analogue of matrix-matrix multiplication.
%
% See also TIMES.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'separableApprox') )           % SEPARABLEAPPROX * ???
    
    if ( isa(g, 'double') )         % SEPARABLEAPPROX * DOUBLE
        if ( numel(g) == 1 )
            h = f;
            if (g ~= 0)
                h.pivotValues = h.pivotValues ./ g;
            else
                h.cols = 0*h.cols(:,1);
                h.rows = 0*h.rows(:,1);
                h.pivotValues = inf;
                if ~isempty( h.pivotLocations )
                    h.pivotLocations = h.pivotLocations(1,:);
                end
            end
        else
            error('CHEBFUN:SEPARABLEAPPROX:mtimes:size', 'Sizes are inconsistent.');
        end
        
    elseif ( isa(g, 'chebfun') )    % SEPARABLEAPPROX * CHEBFUN
        cols = f.cols;
        rows = f.rows;
        fScl = diag( 1./f.pivotValues );
        X = innerProduct( rows, g );
        h = cols * fScl * X;
        
    elseif ( isa(g, 'separableApprox') )   % SEPARABLEAPPROX * SEPARABLEAPPROX compute

        % Type check:
        if ( ~strcmp( class(f),class(g) ) )
            error( 'CHEBFUN:SEPARABLEAPPROX:plus:unknown', ...
                ['Undefined function ''mtimes'' for input arguments of type %s ' ...
                'and %s. Try converting to the same type.'], class(f), class(g));
        end 
                
        % Get the columns and rows of f
        fCols = f.cols;
        fRows = f.rows;
        f_PivotSigns = diag( sign( f.pivotValues ) ); 
        fScl = diag( sqrt( 1./ abs( f.pivotValues ) ) );
        % Share out scaling: 
        fCols = fCols * fScl * f_PivotSigns;
        fRows = fRows * fScl; 
        
        % Get the columns and rows of g
        gCols = g.cols;
        gRows = g.rows;
        g_PivotSigns = diag( sign( g.pivotValues ) ); 
        gScl = diag(sqrt( 1./ abs( g.pivotValues ) ) );
        % Share out scaling: 
        gCols = gCols * gScl;
        gRows = gRows * gScl * g_PivotSigns;
        
        % Compute integral in s.
        X = innerProduct( fRows, gCols );
        
        % Construct low rank form of the result:
        [U, S, V] = svd( X, 'econ' );
        h = f;
        h.cols = fCols * U;
        h.rows = gRows * V;
        h.pivotValues = 1 ./ diag( S );
        
    elseif ( isa(g, 'chebfun2v') )  % SEPARABLEAPPROX * CHEBFUN2V
        nG = g.nComponents; 
        h = g; 
        gc = g.components; 
        for jj = 1:nG 
           h.components{jj} = times(f, gc{jj}); 
        end
        
    else
        error('CHEBFUN:SEPARABLEAPPROX:mtimes:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
        
    end
    
else
    
    h = mtimes(g.', f.').';
    
end

end
