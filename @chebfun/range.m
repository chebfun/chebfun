function r = range(f,dim)
% RANGE   Range of chebfun.
%
% R = RANGE(F) returns the range R = max(F)-min(F) of the chebfun F.
%
% R = RANGE(F,DIM) operates along the dimension DIM of the quasimatrix F. 
% If DIM represents the continuous variable, then R is a vector.
% If DIM represents the discrete dimension, then R is a quasimatrix.
% The default for DIM is 1, unless F has a singleton dimension, in which 
% case DIM is the continuous variable. 
%
% See also CHEBFUN/MINANDMAX.

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ~isreal(f)
    f = abs(f);
end

if nargin == 1
    r = minandmax(f);
else
    r = minandmax(f,dim);
end

if isnumeric(r)
    if numel(f) == 1,
        r = diff(r.');  
        return
    end
    if ~f(1).trans, 
        r = diff(r);    
    else   
        r = diff(r,1,2); 
    end
    return
end

if ~r(1).trans
    r = r(:,2) - r(:,1);
else
    r = r(2,:) - r(1,:);
end

end
