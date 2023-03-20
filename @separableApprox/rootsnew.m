function r = roots( f, g, varargin )
%ROOTS   Zero contours of a SEPARABLEAPPROX.
%   R = ROOTS(F) returns the zero contours of F as a quasimatrix of chebfuns.
%   Each column of R is one zero contour. This command only finds contours when
%   there is a change of sign and it can also group intersecting contours in a
%   non-optimal way.  This command may be inaccurate for components near
%   the boundary of the domain. 
%
%   R = ROOTS(F, G) returns the isolated mutual roots of F and G.
%
%   R = ROOTS(F, G, METHOD) allows the underlying rootfinding algorithm to
%   be supplied. If METHOD = 'ms' or 'marchingsquares', the Marching
%   Squares algorithm is employed, which is fast but not very 
%   robust. If METHOD = 'resultant', a hidden variable resultant method
%   based on Bezout resultants is employed, slower but more robust.
%   See the CHEBFUN2V/ROOTS documentation to see which algorithm is used
%   when no METHOD is passed.
%  
% See also CHEBFUN2V/ROOTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty chebfun input:
if ( isempty( f ) )
    r = []; 
    return
end 

dom = f.domain;

if ( nargin == 1 )  % Seek zero curves of scalar function
    
    if ( length( f ) == 1 )  % Special case: f has rank 1:
        cols = f.cols;
        rows = f.rows;
        yrts = 1i*(roots( cols )+realmin);  
        xrts = roots( rows ) + realmin*1i; % Add complex to ensure it's not real
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
        
    elseif ( isreal( f ) )  % Main case: zero curves for real function 
        % We seek accurate chebfun representations of each component
        % of the zero curve.  First we call Matlab's contourc function
        % to get data points, which typically will lie on the zero curve 
        % to 5-6 digit accuracy and will be smoothly spaced along it 
        % to 2-3 digit accuracy.  The choice n = 502 is used to get a
        % grid that does not involve the boundary of the domain.
        n = 502;
        x = linspace( dom(1), dom(2), n );
        x(1) = []; x(end) = []; 
        y = linspace( dom(3), dom(4), n );
        y(1) = []; y(end) = []; 
        [xx, yy] = meshgrid( x, y );
        vals = feval( f, xx, yy );
        C = contourc( x, y, vals, 0*[1 1] );
        [fx, fy] = grad(f);        
        gradf = @(data) feval( fx, real(data), imag(data)) ...;
                   + 1i*feval( fy, real(data), imag(data));
        fval = @(data) feval( f, real(data), imag(data));
        
        % The contruction of chebfuns proceeds component by component,
        % using complex arithmetic for convenience.
        j = 1; r = chebfun;
        while ( j < length(C) )
            k = j + C(2, j);
            D = C(:, j+1:k);
            data = ( D(1, :) + 1i*(D(2, :)+realmin) ).';            
            for kk = 1:5
                ep = 10^(-3*kk);
                curve = chebfun(data);                   % make chebfun
                curve = simplify(curve,ep);              % simplify it
                data = curve(chebpts(2*length(curve)));  % sample finely
                ff = fval(data);                         % function values
                g = gradf(data);                         % gradient values
                data = data - ff.*( g./abs(g).^2 );      % Newton step
            end
            j = k + 1;
            r = [ r , curve ];
        end

    else  % Function is complex-valued: reduce to real case
        r = roots( [ real(f) ; imag(f) ] );
        if ( ~isempty( r ) )
            r = r(:, 1) + 1i*r(:, 2);
        end
    end
    
elseif ( isa(g, 'separableApprox') )  % seek zero points of vector function
    % TODO: This should not call chebfun2v directly.
    r = roots( chebfun2v( f, g ), varargin{:} );
    
end

end
