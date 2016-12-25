function H = mldivide(C, F)
%\   SPHEREFUNV left divide.
%   C\F Divides each component of the SPHEREFUNV F by the DOUBLE C. 
%
%   Only allowed to divide a SPHEREFUNV by a DOUBLE.
%
% See also SPHEREFUNV/MRDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(C) || isempty(F) )
    H = spherefunv;
    return
end

if ( ~isa(C, 'double') )
    error('SPHEREFUN:SPHEREFUNV:mldivide:nonScalar', ...
        'Division must be by a scalar.');
end

% Left divide.
H = F;
for j = 1 : 3
    H.components{j} = mldivide(C, F.components{j});
end

end
