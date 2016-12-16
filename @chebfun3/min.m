function h = min(f, g, dim)
%MIN   Minimum value of a CHEBFUN3 in one direction.
%   MIN(F) returns a CHEBFUN2 representing the minimum of the CHEBFUN3
%   object F along the x direction, i.e, MIN(F) = @(y,z) min(F(:, y, z)).
%
%   MIN(F, [], DIM) returns a CHEBFUN2 representing the minimum of F along 
%   the dimension DIM. DIM should be 1, 2 or 3 to minimize over X, Y or Z,
%   respectively.
%
%   WARNING: This function is not always accurate to the expected precision.
% 
% See also CHEBFUN3/MIN2 and CHEBFUN3/MIN3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) )
    error('CHEBFUN:CHEBFUN3:min:input', 'CHEBFUN3 is empty');
end

% Default to min of one chebfun3:
if ( nargin < 2 )
    g = [];
end

% Default to minimum along the x direction: 
if ( nargin < 3 )
    dim = 1;
end

% Do not allow min(F, G):
if ( nargin > 1 && ~isempty(g) )
    error('CHEBFUN:CHEBFUN3:min:twoCHEBFUN3Inputs', ...
        'Unable to minimize two CHEBFUN3 objects.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We have no idea how to achieve this in an efficient way. This
% is an attempt to return a result, but typically it won't be accurate to 
% more than 4-5 digits. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dom = f.domain;
n = 129;
if ( dim == 1 )
    vals = sample(f, n, n, n); 
    h = chebfun2(squeeze(min(vals, [], 1)), dom(3:6));
    h = simplify(h);
    
elseif ( dim == 2 )
    vals = sample(f, n, n, n);
    h = chebfun2(squeeze(min(vals, [], 2)), [dom(1:2), dom(5:6)]);
    h = simplify(h);
    
elseif ( dim == 3 )
    vals = sample(f, n, n, n);  
    h = chebfun2(min(vals, [], 3), dom(1:4));
    h = simplify(h);
    
elseif ( dim == 0 )
    error('CHEBFUN:CHEBFUN3:min:dim', 'dim must be 1, 2, or 3.')
else
   % Return the CHEBFUN3 object. This is analogous to that MIN() command in
   % MATLAB.
   h = f;  
end

end