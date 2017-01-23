function f = plus(f, g)
%+   Addition of two CLASSICFUN objects.
%   F + G adds F and G, where F and G may be CLASSICFUN objects or scalars.
%
%   If F and G are both CLASSICFUN objects, they are assumed to have the same domain,
%   and the same mapping. The method gives no warning if their domains/mappings
%   don't agree, but the output of the method will not be useful.
%
% See also MINUS, UPLUS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % CLASSICFUN + [] = []
    
    f = [];

elseif ( isa(g, 'double') ) % CLASSICFUN + double
    
    % Update the onefun, everything else stays the same. The onefun PLUS method
    % ensures that the dimensions of the CLASSICFUN argument and the double argument
    % match (that is, their number of columns):
    f.onefun = f.onefun + g; 
    
elseif ( isa(f, 'double') ) % double + CLASSICFUN
    
    % Call the ONEFUN/PLUS() method for g.onefun:
    g.onefun = g.onefun + f;
    
    % Return g, by assigning it to f:
    f = g;
    
else % CLASSICFUN + CLASSICFUN
    
    % Domains and mappings are assumed to match. Hence, we just need to add the
    % ONEFUNs:
    f.onefun = f.onefun + g.onefun;
      
end

end

