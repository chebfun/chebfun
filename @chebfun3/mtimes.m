function h = mtimes(f, g, varargin)
%*	   Pointwise multiplication for CHEBFUN3 objects.
%
%   c*F or F*c multiplies a CHEBFUN3 F by a scalar c.
%
% See also TIMES.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'chebfun3') )           % CHEBFUN3 * ???
    
    if ( isa(g, 'double') )         % CHEBFUN3 * DOUBLE
        if ( numel(g) == 1 )
            if ( g == 0 ) % This makes factor quasimatrices also zero and 
                % therefore the rank of h will look better than being just 
                % the same as rank(f).
                h = chebfun3(@(x,y,z) 0, f.domain);
            else
                h = f;
                h.core = h.core .* g;
            end
        else
            error('CHEBFUN:CHEBFUN3:mtimes:size', 'Sizes are inconsistent.');
        end
        
    elseif ( isa(g, 'chebfun') )    % CHEBFUN3 * CHEBFUN: mode 1 for now (like CHEBFUN2 * CHEBFUN ).
        % The output is a CHEBFUN2. This is the continuous analogue of
        % tensor vector multiplication which makes a matrix.
        [fCore, fCols, fRows, fTubes] = tucker(f);
        
        X = innerProduct(fCols, g);
        temp = squeeze(chebfun3.txm(fCore, X.'));
        [U, S, V] = svd(temp);
        h = chebfun2();
        h.cols = fRows * U;
        h.rows = fTubes * V;
        h.pivotValues = 1 ./ diag(S);
        
    elseif ( isa(g, 'chebfun2') )   % CHEBFUN3 * CHEBFUN2: modes 1 for now.
        % The output is another CHEBFUN3. Recall that, mode-1 contraction 
        % of an M*N*P tensor by a Q*M matrix creates a Q*N*P tensor. This 
        % is the continuous analogue of that operation.
        [fCore, fCols, fRows, fTubes] = tucker(f);
        
        gRows = g.rows; % rows correspond to the 1st variable of g. This is Chebfun2's convention.
        gCols = g.cols; % cols correspond to the 2nd variable of g. This is Chebfun2's convention.
        gPivots = diag(1./g.pivotValues);
        
        X = innerProduct(fCols, gRows);
        temp = squeeze(chebfun3.txm(fCore, X.'));
        
        % Form the output:
        h = chebfun3();
        h.cols = gCols;   % 2nd variable of g sits in place of 1st variable of h.
        h.rows = fRows;   % 2nd variable of h is the 2nd variable of f.
        h.tubes = fTubes; % 3rd variable of h is 3rd variable of f.
        h.core = chebfun3.txm(temp, gPivots, 1);
        h.domain = f.domain;
                           
    elseif ( isa(g, 'chebfun') )    % CHEBFUN3 * CHEBFUN3
        error('CHEBFUN:CHEBFUN3:mtimes', 'Not written yet.');

    elseif ( isa(g, 'chebfun3v') )  % CHEBFUN3 * CHEBFUN3V
        nG = g.nComponents; 
        h = g; 
        gc = g.components; 
        for jj = 1:nG 
           h.components{jj} = times(f, gc{jj}); 
        end
        
    else
        error('CHEBFUN:CHEBFUN3:mtimes:unknown', ...
            ['Undefined function ''mtimes'' for input arguments of type %s ' ...
            'and %s.'], class(f), class(g));
        
    end
    
elseif ( isa(f, 'double') && isa(g, 'chebfun3') )  % DOUBLE * CHEBFUN3
    h = mtimes(g.', f.').';
else
    error('CHEBFUN:CHEBFUN3:mtimes:size', 'Sizes are inconsistent.');
end

end