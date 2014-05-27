function f=times(f,g)
% .*   Chebfun2 multiplication.
% 
% F.*G multiplies chebfun2 objects F and G.  Alternatively F or G could be
% a double.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if isempty(f) || isempty(g) 
   f = chebfun2;   % just return an empty chebfun2. 
   return; 
end

% ADCHEBFUN2*DOUBLE or DOUBLE*ADCHEBFUN2. Just call mtimes.
if ( isa(f,'double') || isa(g,'double') )

    f = mtimes(f,g);

% ADCHEBFUN2*CHEBFUN2 or CHEBFUN2*ADCHEBFUN2 does not need complicated
% derivative computation, as a CHEBFUN2 does not have any derivative
% information, and can just be treated like a constant
elseif ( isa(f,'ADchebfun2') && isa(g,'chebfun2') ) 

    f.chebfun2 = f.chebfun2.*g;
    
    f.der = f.der*g;

elseif ( isa(f,'chebfun2') && isa(g,'ADchebfun2') )
    
    g.chebfun2 = f.*g.chebfun2;
    
    g.der = (g.der)*f;

    f = g;
elseif isa(f,'chebfun2') && isa(g,'chebfun2v')
%% chebfun2 * chebfun2v

    % This functionality may be taken out of a release.  
    g.xcheb = f.*g.xcheb; 
    g.ycheb = f.*g.ycheb;
    if ~isempty(g.zcheb)
        g.zcheb = f.*g.zcheb;
    end
    f = g;
else

    % We had a chebfun2.*unknown, so complain. 
    error('Chebfun2:times','Can only do chebfun2 times scalar or chebfun2.');
end
end