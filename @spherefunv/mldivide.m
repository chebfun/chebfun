function H = mldivide( f, G )
%\  SPHEREFUNV left divide.
%
%  c\F Divides each component of a SPHEREFUNV by the scalar c. 
%
%  Only allowed to divide a SPHEREFUNV by a scalar.
%
% See also MRDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(f) ) || ( isempty(G) ) )
    H = spherefunv;
    return
end

if ( ~isa(f, 'double') )
    error('SPHEREFUN:SPHEREFUNV:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide.
H = G;
for j = 1 : 3
    H.components{j} = mldivide( f, G.components{j} );
end

end
