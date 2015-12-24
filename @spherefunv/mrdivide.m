function H = mrdivide(F,g)
%/   SPHEREFUNV right divide.
%
% F/c divides each component of a SPHEREFUNV by a scalar. 
% 
% Only allowed to divide by scalars. 
% 
% See also MLDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(F) ) || ( isempty(g) ) )
   H = spherefunv;
   return 
end

if ( ~isa(g,'double') && ~isa(g,'chebfun2') )
    error('SPHEREFUN:SPHEREFUNV:mrdivide:nonScalar', ...
        'Division must be scalar valued.');
end


% componentwise divide. 
if ( isa(g,'double') )
    H = F; 
    for j = 1 : 3
        H.components{j} = mrdivide(F.components{j}, g);
    end
else
    H = F;
    for j = 1 : 3
        H.components{j} = rdivide(F.components{j}, g);
    end
end
    
end
