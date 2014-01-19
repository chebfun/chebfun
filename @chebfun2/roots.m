function r = roots( f, varargin )
%ROOTS zero contours of a chebfun2
%
% R = ROOTS(F), returns the zero contours of F as a quasimatrix of chebfuns.
% Each column of R is one zero contour. This command only finds contours when
% there is a change of sign and it can also group intersecting contours in
% a non-optimal way. Contours are computed to, roughly, four digits of
% precision. In particular, this command cannot reliably compute isolated
% real roots of F.
%
% In the special case when F is of length 1 then the zero contours are
% found to full precision.
%
% R = roots(F,G) returns the isolated points of F and G.
%
% See also CHEBFUN2V/ROOTS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.


accuracy = 1e-5; % Go for five digits.
dom = f.domain;

if ( nargin == 1 ) 
if ( length( f ) == 1 )  % If the Chebfun2 is rank 1:
    cols = f.cols;
    rows = f.rows;
    yrts = roots( cols ) + realmin*1i;  % add complex to ensure its not real
    xrts = roots( rows ) + realmin*1i;
    ry = chebfun( yrts.' );
    rx = chebfun( xrts.' );
    r = [ry rx];
elseif ( isreal( f ) ) 
    % Use Matlab's contourc function.
    n = 2049; % disc size.
    x = linspace( dom(1), dom(2), n ); 
    y = linspace( dom(3), dom(4), n );
    [xx, yy] = meshgrid( x, y ); 
    vals = feval( f, xx, yy );
    C = contourc( x, y, vals, 0*[1 1] );  % Marching squares.
    
    % Store solution in complex-valued chebfuns: 
    j = 1; r = chebfun;
    while ( j < length(C) )
        k = j + C(2, j);  
        D = C(:, j+1:k);
        f = chebfun(D(1,:) +  1i*(D(2,:)+realmin));
        f = simplify( f, accuracy );
        j = k + 1; 
        r = [ r , f ];
    end
else                         % function is complex-valued.
    % Call chebfun2v/roots: 
    r = roots( [ real(f) ; imag(f) ] );  
    % Make them complex again.
    if ( ~isempty( r ) )
        r = r(:, 1) + 1i*r(:, 2);
    end
end
elseif ( isa(varargin{1},'chebfun2') )
    % Bivariate rootfinding: 
    r = roots( chebfun2v( f, varargin{1} ) );
end
    
end