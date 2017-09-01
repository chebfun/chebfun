function r = roots( f, g, varargin )
%ROOTS   Zeros of a SEPARABLEAPPROX.
%   R = ROOTS(F) returns the zero contours of F as a quasimatrix of chebfuns.
%   Each column of R is one zero contour.  This command only finds contours when
%   there is a change of sign, and it can also group intersecting contours in a
%   non-optimal way.
%
%   For a faster plot to graphical accuracy use CONTOUR(F, [0 0]).
%
%   R = ROOTS(F, G) returns the isolated mutual roots of F and G.
%
%   R = ROOTS(F, G, METHOD) allows the underlying rootfinding algorithm to
%   be specified.  If METHOD = 'ms' or 'marchingsquares', the Marching
%   Squares algorithm is employed, which is fast but not very robust.
%   If METHOD = 'resultant', a hidden variable resultant method
%   based on Bezout resultants is employed, slower but more robust.
%   See the CHEBFUN2V/ROOTS documentation to see which algorithm is used
%   when no METHOD is passed.
%  
% Example:
%   cheb.xy;
%   f = x.^2 + y.^2 - 1/4;
%   roots(x,f)                               % [0 -.5; 0 .5]
%   c = roots(f);
%   arclength = sum(abs(diff(c)))            % pi
%   area = abs(sum(real(c).*diff(imag(c))))  % pi/4
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
scl = norm(dom,inf);

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
        % to get initial data points, which typically will lie on the zero
        % curve to 5-6 digit accuracy and will be smoothly spaced along 
        % it to 2-3 digit accuracy.  The choice n = 502 is used to get a
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
        gradf = @(data) feval( fx, real(data), imag(data)) ...
                   + 1i*feval( fy, real(data), imag(data));
        fval = @(data) feval( f, real(data), imag(data));
        
        % The construction of chebfuns proceeds by improving the accuracy
        % of the curves component by component, using complex arithmetic
        % for convenience.
        j = 1; r = chebfun;
        while ( j < length(C) )
            k = j + C(2, j);
            D = C(:, j+1:k);
            data = ( D(1, :) + 1i*(D(2, :)+realmin) ).';            
            ii = find(abs(diff(data)) < 1e-8*scl);
            data(ii) = [];                           % eliminate repetitions
            npts = length(data);
            err = 999; errnew = 1e-2;
            stepno = 0;
            curvenew = chebfun(data);
            while (errnew < err) && (stepno < 6)
                stepno = stepno+1;
                err = errnew;
                curve = curvenew;
                s = [0; cumsum(abs(diff(data)))];    % empirical arclength
                s = 2*s/s(end) - 1;                  % normalize to [-1,1]
                data = interp1(s, data, chebpts ...  % interpolate
                    (length(data)), 'spline');
                curvenew = chebfun(data);            % make chebfun
                data = curvenew(chebpts(npts));      % sample finely
                data = snap(data);                   % snap to boundary
                ff = fval(data);                     % function values
                g = gradf(data);                     % gradient values
                errnew = norm(ff,inf)/vscale(f);     % max rel error
                data = data - ff.*( g./abs(g).^2 );  % Newton step
                curvenew = chebfun(data);            % make chebfun
                curvenew = simplify(curvenew);       % simplify it
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

function d = snap(d)   % adjust endpoints to snap to boundary of domain
    n = length(d);

    if abs(real(d(1))-dom(2)) < .02*scl                 % snap to right bndry
        dx0 = dom(2)-real(d(2)); dx = real(d(1)-d(2));
        d(1) = d(2) + (dx0/dx)*(d(1)-d(2));
    end
    if abs(real(d(n))-dom(2)) < .02*scl
        dx0 = dom(2)-real(d(n-1)); dx = real(d(n)-d(n-1));
        d(n) = d(n-1) + (dx0/dx)*(d(n)-d(n-1));
    end

    if abs(real(d(1))-dom(1)) < .02*scl                 % snap to left bndry
        dx0 = dom(1)-real(d(2)); dx = real(d(1)-d(2));
        d(1) = d(2) + (dx0/dx)*(d(1)-d(2));
    end
    if abs(real(d(n))-dom(1)) < .02*scl
        dx0 = dom(1)-real(d(n-1)); dx = real(d(n)-d(n-1));
        d(n) = d(n-1) + (dx0/dx)*(d(n)-d(n-1));
    end

    if abs(imag(d(1))-dom(4)) < .02*scl                 % snap to top bndry
        dy0 = dom(4)-imag(d(2)); dy = imag(d(1)-d(2));
        d(1) = d(2) + (dy0/dy)*(d(1)-d(2));
    end
    if abs(imag(d(n))-dom(4)) < .02*scl
        dy0 = dom(4)-imag(d(n-1)); dy = imag(d(n)-d(n-1));
        d(n) = d(n-1) + (dy0/dy)*(d(n)-d(n-1));
    end

    if abs(imag(d(1))-dom(3)) < .02*scl                 % snap to bottom bndry
        dy0 = dom(3)-imag(d(2)); dy = imag(d(1)-d(2));
        d(1) = d(2) + (dy0/dy)*(d(1)-d(2));
    end
    if abs(imag(d(n))-dom(3)) < .02*scl
        dy0 = dom(3)-imag(d(n-1)); dy = imag(d(n)-d(n-1));
        d(n) = d(n-1) + (dy0/dy)*(d(n)-d(n-1));
    end

end

end
