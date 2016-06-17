function F = times(F, G)
%.*   Pointwise multiplication of two CHEBFUN3V objects.
%   F.*G computes pointwise, componentwise multiplication of two
%   CHEBFUN3V objects, a CHEBFUN3V and a CHEBFUN3, or a CHBFUN3V 
%   and a double.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) || isempty(G) )
    F = chebfun3v();
    return
end

if ( ~isa(F, 'chebfun3v') ) 
    F = times(G, F); 
    return
end

nF = F.nComponents; 

if ( isa(G, 'double') )             % CHEBFUN3V .* double
    if ( numel(G) == 1 )            % CHEBFUN3V .* scalar
        scalar = G;
        for jj = 1:nF 
            F.components{jj} = times(F.components{jj}, scalar); 
        end
    elseif ( (size(G, 1) == nF) || ( F.isTransposed && (size(G, 2) == nF) ) ) 
                                    % CHEBFUN3V .* vector
        for jj = 1 : nF 
            F.components{jj} = times(F.components{jj}, G(jj));
        end
    else
        error('CHEBFUN:CHEBFUN3V:times:double', ...
            'CHEBFUN3V and double size mismatch.');
    end  
    
elseif ( isa(G, 'chebfun3v') )      % CHEBFUN3V .* CHEBFUN3V
    nG = G.nComponents; 
    if ( nF ~= nG ) 
         error('CHEBFUN:CHEBFUN3V:times:numel', ...
             'CHEBFUN3V components mismatch.');
    end
    for jj = 1:nF 
        F.components{jj} = times(F.components{jj}, G.components{jj}); 
    end
    
elseif ( isa(G, 'chebfun3') )       % CHEBFUN3V .* CHEBFUN3
    for jj = 1 : nF 
        F.components{jj} = times(F.components{jj}, G);
    end
    
else  % error
    error('CHEBFUN:CHEBFUN3V:times:inputs', 'Unrecognized input arguments.');
end

end