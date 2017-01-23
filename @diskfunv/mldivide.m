function G = mldivide( a, F )
%\  DISKFUNV left divide.
%
%  a\F Divides each component of a DISKFUNV F by the scalar a. 
%
%  One is only allowed to divide a DISKFUNV by a scalar.
%
% See also MRDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( ( isempty(a) ) || ( isempty(F) ) )
    G = diskfunv();
    return
end

% Only allow a\F, where a is a scalar:
if ( ~isa(a, 'double') )
    error('DISKFUN:DISKFUNV:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide:
G = F;
G.components{1} = mldivide( a, F.components{1} );
G.components{2} = mldivide( a, F.components{2} );

end