function F = plus( F, G ) 
% + PLUS of two DISKFUNV objects. 
%   F + G if F and G are DISKFUNV objects does componentwise addition. 
%
%   F + G if F is a double and G is a DISKFUNV does componentwise addition. 
% 
%   F + G if F is a DISKFUNV and G is a double does componentwise addition.
% 
%   PLUS(F,G) is called for the syntax F + G. 
% 
% See also DISKFUNV/MINUS 

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) || isempty( G ) ) 
    F = diskfunv; 
    return
end

% Switch order of inputs if F is not a DISKFUNV:
if ( ~isa( F, 'diskfunv' ) )
    F = plus(G, F); 
    return
end

% How many components?:
nF = F.nComponents;

if ( isa(G, 'double') )              % DISKFUNV + DOUBLE
    if ( numel(G) == 1 )             % DISKFUNV + SCALAR
       
       F.components{1} = plus(F.components{1}, G);
       F.components{2} = plus(F.components{2}, G); 
    
    elseif ( numel(G) == nF )        % DISKFUNV + MATRIX
    
        F.components{1} = plus(F.components{1}, G(1));
        F.components{2} = plus(F.components{2}, G(2));         
    
    else
        
        error('CHEBFUN:DISKFUNV:plus:doubleSize', 'Dimension mismatch.')
    
    end
    
elseif ( isa(G, 'diskfun') )        % DISKFUNV + DISKFUN
    
    F.components{1} = plus(F.components{1}, G);
    F.components{2} = plus(F.components{2}, G);  

elseif ( isa(G, 'diskfunv') )       % DISKFUNV + DISKFUNV
    
    if ( G.isTransposed ~= F.isTransposed )
        error('CHEBFUN:DISKFUNV:plus:transposed', 'Dimension mismatch.')
    end
    % Add each component together:
    F.components{1} = plus(F.components{1}, G.components{1});
    F.components{2} = plus(F.components{2}, G.components{2});

else
    
    error('CHEBFUN:DISKFUNV:plus:type', 'Unrecongized input arguments')

end
end