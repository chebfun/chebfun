function H = mldivide( f, G )
%\  DISKFUNV left divide.
%
%  f\F Divides each component of a DISKFUNV by the scalar f. 
%
%  Only allowed to divide a DISKFUNV by a scalar.
%
% See also MRDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
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
H.components{1} = mldivide( f, G.components{1} );
H.components{2} = mldivide( f, G.components{2} );

end
