function h = rdivide(f, g)
%./   Pointwise right divide for CHEBFUN3 objects.
%   RDIVIDE(F, G) is the same as F./G.
%
%   F./G returns (1/G)*F if F is a CHEBFUN3 and G is a double.
%
%   F./G returns F/G if F is a double and G is a CHEBFUN3. This does not 
%   work if G becomes numerically close to zero.
%
% See also CHEBFUN3/LDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) || isempty(g) ) 
    h = chebfun3();
    return
end

if ( isa(f, 'chebfun3') && isa(g, 'chebfun3') )    % CHEBFUN3 ./ CHEBFUN3
    if ( ~domainCheck(f, g))
       error('CHEBFUN:CHEBFUN3:rdivide:domains', 'Domains inconsistent.') 
    end
    h = chebfun3(@(x,y,z) feval(f, x, y, z) ./ feval(g, x, y, z), ...
        f.domain);
    
elseif ( isa(f, 'chebfun3') && isa(g, 'double') )  % CHEBFUN3 ./ double 
    if ( g == 0 )
        error('CHEBFUN:CHEBFUN3:rdivide:divByZero', ...
            'Division by zero or near zero.')
    end
    h = f .* (1/g);
        
elseif ( isa(f, 'double') && isa(g, 'chebfun3') )   
       [ss, zeroVal] = singleSignTest(g);
       if ( ( ss == 1 ) && ( zeroVal == 0 ) )
           h = chebfun3(@(x,y,z) f ./ feval(g, x, y, z), g.domain);
       else
          error('CHEBFUN:CHEBFUN3:rdivide:zero', ...
              'Attempting to invert a CHEBFUN3 with a root.'); 
       end
       
elseif ( isa(f, 'chebfun3') && isa(g, 'chebfun3v') )
    % TODO: RDIVIDE on the components: 
    
else
    error('CHEBFUN:CHEBFUN3:rdivide:badInputs', 'Unrecognised operation.');
    
end

end