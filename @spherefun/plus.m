function h = plus(f, g)
%+   Plus for SPHEREFUN objects.
%
% F + G adds F and G. F and G can be scalars or SPHEREFUN objects.

if ( ~isa(f, 'spherefun') ) % ??? + SPHEREFUN
    
    h = plus(g, f);
    
elseif ( isempty(g) ) % SPHEREFUN + []
    
    h = g; 
    
elseif ( isempty(f) ) % [] + SPHEREFUN
    
    h = f; 
    
elseif ( isa( g, 'double' ) )           % SPHEREFUN + DOUBLE
    
    g = compose( 0*f,@plus, g);   % promote double to object class.  
    h = plus(f, g); 
    
elseif ( ~isa(g, 'spherefun') )          % SPHEREFUN + ???
    
    error( 'SPHEREFUN:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class( f ), class( g ));
        
else                                     % SPHEREFUN + SPHEREFUN
    
    % Domain Check:
    if ( ~domainCheck(f, g) )
        error('SPHEREFUN:plus:domain', 'Inconsistent domains.');
    end
    
    % Check for zero SPHEREFUN objects:
    if ( iszero(f) )
        h = g;
    elseif ( iszero(g) )
        h = f;
    else
        % Add together two nonzero SPHEREFUN objects:
        % The algorithm is as follows: Split f and g into their plus/minus
        % components.  Do the compression_plus algorithm described in
        % @separableApprox/compression plus on each pair of plus and minus
        % components.
        
        [fp,fm] = partition(f);
        [gp,gm] = partition(g);
        
        hp = plus@separableApprox(fp,gp);
        hm = plus@separableApprox(fm,gm);
        
        % Put pieces back together.
        h = combine(hp,hm);
    end 
    
end

end

