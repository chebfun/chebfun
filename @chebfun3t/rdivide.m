function h = rdivide(f, g)
%./   Pointwise right divide for CHEBFUN3T objects.
%   RDIVIDE(F, G) is the same as F./G.
%
%   F./G returns (1/G)*F if F is a CHEBFUN3T and G is a double.
%
%   F./G returns F/G if F is a double and G is a CHEBFUN3T. This does not 
%   work if G becomes numerically close to zero.
%
% See also CHEBFUN3T/LDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) || isempty(g) ) 
    h = chebfun3t();
    return
end

if ( isa(f, 'chebfun3t') && isa(g, 'chebfun3t') ) % CHEBFUN3T ./ CHEBFUN3T
    if ( ~domainCheck(f, g))
       error('CHEBFUN:CHEBFUN3T:rdivide:domains', 'Domains inconsistent.') 
    end
    h = chebfun3t(@(x,y,z) feval(f, x, y, z) ./ feval(g, x, y, z), ...
        f.domain);
    
elseif ( isa(f, 'chebfun3t') && isa(g, 'double') )  % CHEBFUN3T ./ double 
    if ( g == 0 )
        error('CHEBFUN:CHEBFUN3T:rdivide:divByZero', ...
            'Division by zero or near zero.')
    end
    h = f .* (1/g);
        
elseif ( isa(f, 'double') && isa(g, 'chebfun3t') )   
           h = chebfun3t(@(x,y,z) f ./ feval(g, x, y, z), g.domain);
           
else
    error('CHEBFUN:CHEBFUN3T:rdivide:badInputs', 'Unrecognised operation.');
    
end

end