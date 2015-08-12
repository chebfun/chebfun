function h = plus(f, g)
%+   Plus for SPHEREFUN/DISKFUN objects.
%
% F + G adds F and G. F and G can be scalars or DISKFUN objects.

if ( ~isa(f, 'diskfun') ) % ??? + DISKFUN
    
    h = plus(g, f);
    
elseif ( isempty(g) ) % DISKFUN + []
    
    h = g; 
    
elseif ( isempty(f) ) % [] + DISKFUN
    
    h = f; 
    
elseif ( isa( g, 'double' ) )           % DISKFUN + DOUBLE
    
    g = compose( 0*f,@plus, g);   % promote double to object class.  
    h = plus(f, g); 
    
elseif ( ~isa(g, 'diskfun') )          % DISKFUN + ???
    
    error( 'DISKFUN:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class( f ), class( g ));
        
else                                     % DISKFUN + DISKFUN
    
    % Domain Check:
    if ( ~domainCheck(f, g) )
        error('DISKFUN:plus:domain', 'Inconsistent domains.');
    end
    
    % Check for zero DISKFUN objects:
    if ( iszero(f) )
        h = g;
    elseif ( iszero(g) )
        h = f;
    else
        % Add together two nonzero DISKFUN objects:
        % The algorithm is as follows: Split f and g into their plus/minus
        % components.  Do the compression_plus algorithm described in
        % @separableApprox/compression plus on each pair of plus and minus
        % components.
        
        [fp,fm] = partition(f);
        [gp,gm] = partition(g);
        
        hp = plus@separableApprox(fp,gp);
        hp.idxPlus = 1:size(hp.cols,2);
        
        hm = plus@separableApprox(fm,gm);
        hm.idxMinus = 1:size(hm.cols,2);
        
        % Put pieces back together.
        h = combine(hp,hm);
    end 
    
end

end

