function G = mrdivide(F, a)
%/   DISKFUNV right divide.
%
% F/a divides each component of a DISKFUNV F by the scalar a. 
% 
% Only allowed to divide by scalars. 
% 
% See also MLDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( ( isempty(F) ) || ( isempty(a) ) )
   G = diskfunv;
   return 
end

% Only allow F/a, where a is a scalar: 
if ( ~isa(a, 'double') )
    error('CHEBFUN:DISKFUNV:mrdivide:nonScalar', ...
        'Division must be scalar valued.');
end

% Componentwise divide. 
G = F; 
G.components{1} = mrdivide(F.components{1}, a);
G.components{2} = mrdivide(F.components{2}, a);

end