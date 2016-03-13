function H = mldivide( f, G )
%\  DISKFUNV left divide.
%
%  c\F Divides each component of a DISKFUNV by the scalar c. 
%
%  Only allowed to divide a DISKFUNV by a scalar.
%
% See also MRDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(f) ) || ( isempty(G) ) )
    H = diskfunv;
    return
end

if ( ~isa(f, 'double') )
    error('DISKFUN:DISKFUNV:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide.
H = G;
for j = 1:G.nComponents
    H.components{j} = mldivide( f, G.components{j} ); 
end

end
