function H = mldivide( f, G )
%\  CHEBFUN2V left divide.
%
%  c\F Divides each component of a CHEBFUN2V by the scalar c. 
%
%  Only allowed to divide a CHEBFUN2V by a scalar.
%
% See also MRDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(f) ) || ( isempty(G) ) )
    H = chebfun2v;
    return
end

if ( ~isa(f, 'double') )
    error('CHEBFUN:CHEBFUN2V:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide.
H = G;
for j = 1:G.nComponents
    H.components{j} = mldivide( f, G.components{j} ); 
end

end
