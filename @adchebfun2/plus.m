function h = plus ( f, g )
%+	  Plus for ADCHEBFUN2
%
% F + G adds ADCHEBFUN2 objects F and G, or a scalar to a 
% ADCHEBFUN2 if either F or G is a scalar.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty(f) || isempty(g) )  % check for empty ADchebfun2
   h = ADchebfun2;  % just return an empty ADchebfun2. 
   return; 
end

if ( ~isa(f, 'adchebfun2') ) % First argument was not an ADchebfun2
    % Swap arguments.
    h = plus(g, f);
    return
end

if isa(g, 'adchebfun2')     % ADCHEBFUN2 + ADCHEBFUN2
    h = f;
    
    % Add the chebfun2s
    h.chebfun2 = f.chebfun2 + g.chebfun2;
    
    % Add the der fields
    h.der = f.der + g.der;
    
elseif ( isa(g, 'chebfun2') || isa(g, 'double') ) % ADCHEBFUN2 + DOUBLE/SCALAR
    f.chebfun2 = f.chebfun2 + g;
else
    error('ADCHEBFUN2:plus:type','Cannot add these two objects together');
end

end