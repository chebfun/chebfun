function f = plus(f, g)
%+	Addition of two BNDFUN objects.
%   F + G adds F and G, where F and G may be FUNCHEB objects or scalars.
%
% See also MINUS, UPLUS.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(g) ) % BNDFUN + [] = []
    
    f = [];

elseif ( isa(g, 'double') ) % BNDFUN + double
    
    % Update the onefun, everything else stays the same. The onefun plus method
    % ensures that the dimensions of the BNDFUN argument and the double argument
    % match (that is, their number of columns):
    f.onefun = f.onefun + g; 
    
elseif ( isa(f,'double') ) % double + FUNCHEB
    
    % Call the ONEFUN/plus method for g.onefun:
    g.onefun = g.onefun + f;
    
    % Return g, by assigning it to f:
    f = g;
    
else % BNDFUN + BNDFUN
    
    % Make sure both BNDFUN objects have the same domain.
    if ( ~checkDomain(f, g))
        
        % Throw an error if domains don't match:
        error('BNDFUN:plus:domainMismatch',...
            'Bndfuns domains must agree.')
    end

    % Domains match, and since we're only working with BNDFUNs, their linear
    % mappings must match as well. Hence, we just need to add the onefuns:
    f.onefun = f.onefun + g.onefun;
      
end

end

