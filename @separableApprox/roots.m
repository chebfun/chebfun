function r = roots( f, g, varargin )
%ROOTS   Zero contours of a SEPARABLEAPPROX.
%   R = ROOTS(F), returns the zero contours of F as a quasimatrix of chebfuns.
%   Each column of R is one zero contour. This command only finds contours when
%   there is a change of sign and it can also group intersecting contours in a
%   non-optimal way. Contours are computed to, roughly, four digits of
%   precision. In particular, this command cannot reliably compute isolated real
%   roots of F or zero curves lying close to the boundary of the domain. 
%
%   In the special case when F is of length 1 then the zero contours are found
%   to full precision.
%
%   R = ROOTS(F, G) returns the isolated points of F and G.
%
%   R = ROOTS(F, G, METHOD) the underlying rootfinding algorithm can be
%   supplied. If METHOD = 'ms' or METHOD = 'marchingsquares', then the Marching
%   Squares algorithm is employed. The Marching Squares algorithm is fast but
%   not particularly robust. If METHOD = 'resultant', then a hidden variable
%   resultant method based on Bezout resultants is employed. The Resultant
%   method is slower but far more robust. See the CHEBFUN2V/ROOTS()
%   documentation to see which algorithm is used when no METHOD is passed.
%  
% See also CHEBFUN2V/ROOTS.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% Check for empty:
if ( isempty( f ) )
    r = []; 
    return
end 

tol = 1e-5; % Go for five digits.
dom = f.domain;

if ( nargin == 1 )
    
    if ( length( f ) == 1 )  
        % The SEPARABLEAPPROX is rank 1:
        
        cols = f.cols;
        rows = f.rows;
        yrts = 1i*(roots( cols )+realmin);  
        xrts = roots( rows ) + realmin*1i; % Add complex to ensure its not real
        r = chebfun; 
        % Go though col(yrts) = 0, make into a horizontal line: 
        dom = rows.domain;
        len = (dom(2)-dom(1))/2; 
        for j = 1 : numel(yrts)
            f = chebfun(@(x) (len*(x+1)+dom(1)) + yrts(j));
            r = [r f];
        end
        % Go though row(xrts) = 0, make into a vertical line: 
        dom = cols.domain;
        len = (dom(2)-dom(1))/2; 
        for j = 1 : numel(xrts)
            f = chebfun(@(x) 1i*(len*(x+1)+dom(1)) + xrts(j));
            r = [r f];
        end
        
    elseif ( isreal( f ) )
        % Function is real-valued.
        
        % Use Matlab's contourc function (Marching Squares). Note n = 502 is
        % chosen so that use a grid that does not involve the boundary of the
        % domain.
        n = 502; % disc size.
        x = linspace( dom(1), dom(2), n );
        x(1) = []; x(end) = []; 
        y = linspace( dom(3), dom(4), n );
        y(1) = []; y(end) = []; 
        [xx, yy] = meshgrid( x, y );
        vals = feval( f, xx, yy );
        C = contourc( x, y, vals, 0*[1 1] );
        
        % Store solution in complex-valued CHEBFUNs:
        j = 1; r = chebfun;
        while ( j < length(C) )
            k = j + C(2, j);
            D = C(:, j+1:k);
            f = chebfun( (D(1,:) +  1i*(D(2,:)+realmin)).' );
            f = simplify( f, tol, 'globaltol' );
            j = k + 1;
            r = [ r , f ];
        end
    else
        % Function is complex-valued.
        
        % Call CHEBFUN2V/ROOTS():
        r = roots( [ real(f) ; imag(f) ] );
        % Make them complex again.
        if ( ~isempty( r ) )
            r = r(:, 1) + 1i*r(:, 2);
        end
    end
    
elseif ( isa(g, 'separableApprox') )
    % TODO: This should not call chebfun2v directly.
    r = roots( chebfun2v( f, g ), varargin{:} );
    
end

end
