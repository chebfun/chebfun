function F = vertcat(F,G)
%VERTCAT Vertical concatenation of chebfun2v objects.
% 
% [F;f] where F is a chebfun2v with two components, and f is a chebfun2
% or scalar then returns a chebfun2v with three components.  The first
% and second component remain unchanged and the third component is f. 
% 
% [f;F] where F is a chebfun2v with two components, and f is a chebfun2 or
% scalar then returns a chebfun2v with three components. The first is f,
% and the second and third are the first and second components of F. 

% Copyright 2013 by The University of Oxford and The Chebfun2 Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun2 information.

if isa(G,'double') && isempty(F.zcheb)
    rect = F.xcheb.corners; 
    F.zcheb = chebfun2(G,rect); 
    return; 
elseif isa(G,'chebfun2') && isempty(F.zcheb)
    rect = F.xcheb.corners;
    rectcheck = G.corners; 
    if all(rect - rectcheck)
       error('CHEBFUN2V:VERTCAT','Chebfun2 object must be on the same domain.'); 
    end
    F.zcheb = G;
    return; 
elseif isa(F,'double') && isempty(G.zcheb)
    rect = G.xcheb.corners; 
    firstcomponent = chebfun2(G,rect); 
    G.zcheb = G.ycheb; 
    G.ycheb = G.xcheb; 
    G.xcheb = firstcomponent; 
    F = G; return; 
elseif isa(F,'chebfun2') && isempty(G.zcheb)
    rect = G.xcheb.corners;
    rectcheck = F.corners; 
    if all(rect - rectcheck)
       error('CHEBFUN2V:VERTCAT','Chebfun2 object must be on the same domain.'); 
    end
    firstcomponent = chebfun2(F,rect); 
    G.zcheb = G.ycheb; 
    G.ycheb = G.xcheb; 
    G.xcheb = firstcomponent;
    F = G; return; 
else
    error('CHEBFUN2V:VERTCAT','Vertical concatenation of these objects is not supported.')
end

end