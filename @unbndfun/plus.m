function f = plus(f, g)
%+	Addition of two UNBNDFUN objects.
%   F + G adds F and G, where F and G must be unbndfun objects or scalars.
%
% See also MINUS, UPLUS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % UNBNDFUN + [] = []
    
    f = [];

elseif ( isa(g, 'double') ) % UNBNDFUN + double
    
% Call the ONEFUN/plus method for f.onefun:
    f.onefun = f.onefun + g; 
    
elseif ( isa(f,'double') ) % double + UNBNDFUN
    
    % Call the ONEFUN/plus method for g.onefun:
    g.onefun = g.onefun + f;
    
    % Return g, by assigning it to f:
    f = g;
    
else % UNBNDFUN + UNBNDFUN
    
    % Make sure UNBNDFUN objects f and g have the same domain.
    if ( ~checkDomain(f, g))
        
        % Throw an error if domains do not match:
        error('CHEBFUN:UNBNDFUN:plus:domainMismatch',...
            'Unbndfuns domains must agree.')
    end

    % Domains match, and since we're only working with BNDFUNs, their linear
    % mappings must match as well. Hence, we just need to add the onefuns:
    f.onefun = f.onefun + g.onefun;
      
end

end