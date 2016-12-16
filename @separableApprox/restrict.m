function f = restrict(f, dom)
%RESTRICT   Restrict the domain of a SEPARABLEAPPROX.
%
% F = RESTRICT(F, DOM) returns a SEPARABLEAPPROX on the domain DOM that
% approximates F on that domain.  DOM should be a vector of length 4 giving the
% coordinates of the corners.

% Copyright 2017 by The University of Oxford and The Chebfun2 Developers.
% See http://www.chebfun.org/ for Chebfun2 information.

if ( isa( dom, 'double' ) )    
    if ( numel( dom ) == 4 )                   % Restrict to DOM. 
        xlen = diff( dom(1:2) );
        ylen = diff( dom(3:4) );
        
        if ( ( xlen == 0 ) && ( ylen == 0) )   % DOM is a point.
            f = feval(f, dom(1), dom(3));
        elseif ( xlen == 0 )                   % DOM is a vertical line
            cols = restrict(f.cols, dom(3:4));
            rows = feval(f.rows, dom(1)); 
            d = 1./f.pivotValues; 
            % Set infinite values to zero.
            d(d == inf) = 0;                   
            f = cols * diag( d ) * rows.';  
        elseif ( ylen == 0 )                   % DOM is a horizontal line
            rows = restrict(f.rows, dom(1:2));
            cols = feval(f.cols, dom(3)); 
            d = 1./f.pivotValues; 
            % Set infinite values to zero.
            d(d == inf) = 0;  
            f = cols * diag( d ) * rows.'; 
        else                                   % DOM is not degenerate
            f.cols = restrict(f.cols, dom(3:4));
            f.rows = restrict(f.rows, dom(1:2));
            f.domain = dom;
        end
    else
        error('CHEBFUN:SEPARABLEAPPROX:restrict:domain', 'Domain not determined.');
    end
    
elseif ( isa(dom, 'chebfun') )
    f = compose(dom, f);
    
else
    error('CHEBFUN:SEPARABLEAPPROX:restrict:domain', 'Unrecognizable domain.');
end

end
