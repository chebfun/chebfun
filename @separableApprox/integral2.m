function I = integral2( f, c )
%INTEGRAL2  Double integral of a SEPARABLEAPPROX over its domain.
%   I = INTEGRAL2(F) returns a value representing the double integral of a
%   SEPARABLEAPPROX.
%
%   I = INTEGRAL2(F, [a b c d]) integrate F over the rectangle region [a b] x [c
%   d] provide this rectangle is in the domain of F.
%
%   I = INTEGRAL2(F, C) computes the volume under the surface F over the region
%   D with boundary C. C should be a complex-valued CHEBFUN that represents a
%   closed curve. This can be a very slow feature.
%
% See also INTEGRAL, SUM2, QUAD2D.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) ) 
    I = 0;
    return
end

% Get the low rank representation for f. 
[cols, D, rows] = cdr(f);

if ( nargin == 1 )
    
    % Double definite integral:
    I = sum( cols ) * D * sum( rows ).';
    
elseif ( nargin == 2 )
    
    if ( isa( c, 'chebfun' ) ) 
        % Integral in area enclosed by CHEBFUN.
        
        if ( ~isreal( c ) )
            % Use Green's theorem to do integration over the domain contained by
            % the CHEBFUN.
            
            dom = f.domain;
            
            % Green's theorem tells that you can integrate a function over a
            % region by integration along the boundary of the region's boundary.
            % 
            % HOW IT WORKS: Green's theorem states that: 
            % 
            %  (*) int_c dot((M;L), normal(c)) = intt_D div( dM/dx ; dL/dy)dxdy. 
            % 
            % The algorithm below writes f = div( dM/dx ; dL/dy) and constructs
            % functions M(x,y) and L(x,y) so the equality (*) holds.
            %
            % There are plenty of different ways of applying Green's theorem. 
            % AT's choice is to treat both variables equally,
            %
            % op = @(x,y) sum( chebfun(@(s) feval(f, s*x, s*y).*s, [0 1] ) );
            % Fs = chebfun2(op , dom, 'vectorize');
            % x = chebfun2( @(x,y) x, dom ); 
            % y = chebfun2( @(x,y) y, dom );
            % F = [ -Fs.*y; Fs.*x ];
            % I = integral(F, c);
            %
            % Another choice would be L = 0 and M = cumsum(f, 2). 
            %
            % M = cumsum(f, 2); 
            % zeroCheb = chebfun2( @(x,y) 0*x, dom );
            % F = [ M ; zeroCheb ];
            % I = integral(F, c);
            
            % [TODO: Is this a good choice?]
            
            op = @(x,y) sum( chebfun(@(s) feval(f, s*x, s*y).*s, [0 1] ) );
            Fs = chebfun2(op , dom, 'vectorize');
            x = chebfun2( @(x,y) x, dom ); 
            y = chebfun2( @(x,y) y, dom );
            F = [ -Fs.*y; Fs.*x ];
            I = integral(F, c);
            
        else
            error('CHEBFUN:SEPARABLEAPPROX:integral2:input', ...
                'Integration path must be complex-valued');
        end
        
    elseif ( isa( c, 'double' ) )  
        % Integral over restricted rectangle
        
        restriction = c;
        if ( length( restriction ) == 4 )   
            % Calculate integral over restriction rectangle.
            g = restrict(f, restriction);
            if ( isa(g, 'separableApprox') )
                I = integral2( restrict( f, restriction ) );
            elseif ( isa(g, 'chebfun') )
                I = sum( restrict ( f,restriction ) );
            end
        else
            error('CHEBFUN:SEPARABLEAPPROX:integral2:baddomain', ...
            'Domain should have four corners.');
        end
        
    else
        error('CHEBFUN:SEPARABLEAPPROX:integral2:nargin', 'Too many input arguments.');
    end
    
end

end
