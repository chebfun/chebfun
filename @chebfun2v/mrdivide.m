function H = mrdivide(F,g)
%/   CHEBFUN2V right divide.
%
% F/c divides each component of a CHEBFUN2V by a scalar. 
% 
% Only allowed to divide by scalars. 
% 
% See also MLDIVIDE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ( isempty(F) ) || ( isempty(g) ) )
   H = chebfun2v;
   return 
end

if ( ~isa(g,'double') && ~isa(g,'chebfun2') )
    error('CHEBFUN:CHEBFUN2V:mrdivide:nonScalar', ...
        'Division must be scalar valued.');
end


% componentwise divide. 
if ( isa(g,'double') )
    H = F; 
    for j = 1 : F.nComponents
        H.components{j} = mrdivide(F.components{j}, g);
    end
else
    H = F;
    for j = 1 : F.nComponents
        H.components{j} = rdivide(F.components{j}, g);
    end
end
    
end
