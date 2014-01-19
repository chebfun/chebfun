function f = times( f , g ) 
%.* Times of two chebfun2v objects. 
%
% F.*G if F is a chebfun2v and G is double returns the chebfun2v after
% componentwise multiplication. 
%
% F.*G if F is a double and G is a chebfun2v returns the chebfun2v after
% componentwise multiplication.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isa(f,'double') ) 
%% double.*chebfun2v
    if ( numel(f) == 1 ) 
        % Assume [f f].'.*g;
        const = f; f = g; 
        f.xcheb = times(g.xcheb,const); 
        f.ycheb = times(g.ycheb,const); 
        f.zcheb = times(g.zcheb,const); 
    elseif ( numel(f) == 2) && isempty(f.zcheb)
        if ( size(f,2) == 1 ) 
            const = f; f=g; 
            f.xcheb = times(g.xcheb,const(1)); 
            f.ycheb = times(g.ycheb,const(2)); 
        else
            error('CHEBFUN2v:times:double','Chebfun2v and double size mismatch.');
        end
    elseif ( numel(f) == 3) && ~isempty(g.zcheb)
        if ( size(f,2) == 1 ) 
           const = f; f=g; 
           f.xcheb = times(g.xcheb,const(1)); 
           f.ycheb = times(g.ycheb,const(2)); 
           f.zcheb = times(g.zcheb,const(3));         
        else
            error('CHEBFUN2v:times:double','Chebfun2v and double size mismatch.');
        end  
    else
        error('CHEBFUN2v:times:double','Chebfun2v times vector of length more than 2.');
    end
elseif ( isa(g,'double') ) 
%% chebfun2v.*double 
    f = times(g,f);
elseif ( isa(f,'chebfun2v') && isa(g,'chebfun2v') )
%% chebfun2v .* chebfun2v
    % componentwise multiplication. 
    if ( isempty(f.zcheb) && isempty(g.zcheb) ) || ( ~isempty(f.zcheb) && ~isempty(g.zcheb) )
        f.xcheb = f.xcheb .* g.xcheb; 
        f.ycheb = f.ycheb .* g.ycheb;
        f.zcheb = f.zcheb .* g.zcheb;
    else
        error('CHEBFUN2v:times','Chebfun2v size mismatch.');
    end
elseif isa(f,'chebfun2v') && isa(g,'chebfun2')
%% chebfun2 * chebfun2v

    % This functionality may be taken out of a release.  
    f.xcheb = g.*f.xcheb; 
    f.ycheb = g.*f.ycheb;
    if ~isempty(f.zcheb)
        f.zcheb = g.*f.zcheb;
    end
else  % error
    error('CHEBFUN2v:times:inputs','Chebfun2v can only times to a double');
end
end