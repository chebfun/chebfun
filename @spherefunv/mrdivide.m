function H = mrdivide(F, G)
%/   SPHEREFUNV right divide.
%   F/G divides each component of the SPHEREFUNV F by the DOUBLE or 
%   SPHEREFUN G.
% 
%   Only allowed to divide by a DOUBLE or a SPHEREFUN.
% 
% See also SPHEREFUNV/MLDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) || isempty(G) )
   H = spherefunv;
   return 
end

if ( ~isa(G,'double') && ~isa(G,'chebfun2') )
    error('SPHEREFUN:SPHEREFUNV:mrdivide:nonScalar', ...
        'Division must be scalar valued.');
end


% componentwise divide. 
if ( isa(G,'double') )
    H = F; 
    for j = 1:3
        H.components{j} = mrdivide(F.components{j}, G);
    end
else
    H = F;
    for j = 1:3
        H.components{j} = rdivide(F.components{j}, G);
    end
end
    
end
