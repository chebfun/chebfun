function f = plus(f, g)
%+   Addition of two BNDFUN objects.
%   F + G adds F and G, where F and G may be BNDFUN objects or scalars.
%
%   If F and G are both BNDFUN objects, they are assumed to have the same
%   domain. The method gives no warning if their domains don't agree, but the
%   output of the method will be gibberish.
%
% See also MINUS, UPLUS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % BNDFUN + [] = []
    
    f = [];

elseif ( isa(g, 'double') ) % BNDFUN + double
    
    % Update the onefun, everything else stays the same. The onefun PLUS method
    % ensures that the dimensions of the BNDFUN argument and the double argument
    % match (that is, their number of columns):
    f.onefun = f.onefun + g; 
    
elseif ( isa(f,'double') ) % double + BNDFUN
    
    % Call the ONEFUN/plus method for g.onefun:
    g.onefun = g.onefun + f;
    
    % Return g, by assigning it to f:
    f = g;
    
else % BNDFUN + BNDFUN
    
    % Domains are assumed to match, and since we're only working with BNDFUNs,
    % their linear mappings must match as well. Hence, we just need to add the
    % onefuns:
    f.onefun = f.onefun + g.onefun;
      
end

end

