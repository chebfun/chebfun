function H = mldivide( f, G )
%\  Chebfun2v left divide.
%
%  c\F Divides each component of a chebfun2v by the scalar c. 
%
%  Only allowed to divide a chebfun2v by a scalar.
%
% See also MRDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( ( isempty(f) ) || ( isempty(G) ) )
    H = chebfun2v;
    return
end

if ( ~isa(f,'double') )
    error('CHEBFUN2:MRDIVIDE:NONSCALAR','Division must be by a scalar.');
end

% Left divide.
H = G;
for j = 1:G.ncomponents
    H.components{j} = mldivide( f, G.components{j} ); 
end

end