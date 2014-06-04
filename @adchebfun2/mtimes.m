function h = mtimes(f,g)
%*	ADchebfun2 multiplication.
%
% c*F or F*c multiplies an ADchebfun2 F by a scalar c.
%
% See also TIMES.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(f) || isempty(g) )  % just return an empty chebfun2
    h = chebfun2;
    return;
end

if ( isa(f,'double') )
    
        % double * chebfun2
        h = g;
    
        % Multiply the chebfun2
        h.chebfun2 = f * ( vertcat( g.chebfun2  ));
       
        % Multiply the derivative
        h.der = ( g.der ) * f;
    
elseif ( isa(g,'double') )
    %% chebfun2 * double
    h = f;
    
    % Multiply the chebfun2
    h.chebfun2 = g*(f.chebfun2);
    
    % Multiply the derivative
    h.der = (f.der)*g;
end