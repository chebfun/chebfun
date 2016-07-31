function H = mrdivide(F,g)
%/   DISKFUNV right divide.
%
% F/g divides each component of a DISKFUNV by a scalar g. 
% 
% Only allowed to divide by scalars. 
% 
% See also MLDIVIDE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(F) ) || ( isempty(g) ) )
   H = diskfunv;
   return 
end

if ( ~isa(g,'double') )
    error('CHEBFUN:DISKFUNV:mrdivide:nonScalar', ...
        'Division must be scalar valued.');
end
% componentwise divide. 
H.components{1} = mrdivide(F.components{1}, g);
H.components{2} = mrdivide(F.components{2}, g);


end
