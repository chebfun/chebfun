function F = times( F , G ) 
%.*   Componentwise multiplication for DISKFUNV objects. 
%   F.*G if F is a DISKFUNV and G is double returns the DISKFUNV after
%   componentwise multiplication.
%
%   F.*G if F is a double and G is a DISKFUNV returns the DISKFUNV after
%   componentwise multiplication.
% 
% See also DISKFUNV/MTIMES

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) || isempty( G ) )
    F = diskfunv;
    return
end

% Switch order of arguments if the first is not a DISKFUNV:
if ( ~isa(F, 'diskfunv') ) 
    F = times(G, F); 
    return
end

% How many components?:
nF = F.nComponents; 

if ( isa(G, 'double') )             % DISKFUNV.*double
    
    if ( numel(G) == 1 )            % DISKFUNV.*scalar
        scalar = G;
        F.components{1} = times(F.components{1}, scalar); 
        F.components{2} = times(F.components{2}, scalar);
        
    elseif ( (size(G, 1) == nF) || ( F.isTransposed && (size(G, 2) == nF) ) )
        
        F.components{1} = times( F.components{1}, G(1) );
        F.components{2} = times( F.components{1}, G(2) );    
        
    else
        
        error('CHEBFUN:DISKFUNV:times:double', ...
            'DISKFUNV and double size mismatch.');
        
    end  
    
elseif ( isa(G, 'diskfunv') )      % DISKFUNV . * DISKFUNV
    
    F.components{1} = times(F.components{1}, G.components{1}); 
    F.components{2} = times(F.components{2}, G.components{2}); 
    
elseif ( isa(G, 'diskfun') )       % DISKFUN * DISKFUNV
    
     F.components{1} = times(F.components{1}, G); 
     F.components{2} = times(F.components{2}, G); 
    
else  % error
    
    error( 'CHEBFUN:DISKFUNV:times:inputs', 'Unrecognized input arguments.' );
    
end
end