function f = plus( f , g ) 
% + PLUS of two chebfun2v objects.

if ( isa(f, 'chebfun2v') )     % CHEBFUN2V + ??? 
    
    if ( isa(g, 'double') )    % CHEBFUN2V + DOUBLE 
    


if ( isa(f,'double') )         % double +
    if(numel(f) == 1) 
        % Assume g + [f f].';
        const = f; f= g; 
        f.xcheb = plus(g.xcheb,const); 
        f.ycheb = plus(g.ycheb,const); 
    elseif( numel(f) == 2) 
        const = f; f=g; 
        f.xcheb = plus(g.xcheb,const(1)); 
        f.ycheb = plus(g.ycheb,const(2)); 
    else
        error('CHEBFUN2v:plus:double','Chebfun2v plus vector of length more than 2.');
    end
elseif( isa(g,'double') ) % chebfun2v + double 
    if(numel(g) == 1) 
        % Assume f + [g g].';
        const = g; temp = f; 
        temp.xcheb = plus(f.xcheb,const); 
        temp.ycheb = plus(f.ycheb,const); 
        f=temp; 
    elseif( numel(g) == 2 ) 
        const = g; temp=f; 
        temp.xcheb = plus(g.xcheb,const(1)); 
        temp.ycheb = plus(g.ycheb,const(2)); 
        f=temp; 
    else
        error('CHEBFUN2v:plus:double','Chebfun2v plus vector of length more than 2.');
    end
elseif (isa(f,'chebfun2v') && isa(g,'chebfun2v') )  % chebfun2v + chefun2v
    temp = f; 
    % plus componentwise. 
    temp.xcheb = plus(f.xcheb,g.xcheb); 
    temp.ycheb = plus(f.ycheb,g.ycheb); 
    f=temp; 
else  % error
    error('CHEBFUN2v:plus:inputs','Chebfun2v can only plus to chebfun2v or double');
end
end