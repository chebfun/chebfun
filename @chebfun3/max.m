function h = max(f, g, dim)
%MAX   Maximum value of a CHEBFUN3 in one of the directions.
%   MAX(F) returns a CHEBFUN2 representing the maximum of the CHEBFUN3 
%   object F along the x direction, i.e, MAX(F) = @(y,z) max(F(:, y, z)).
%
%   MAX(F, [], dim) returns a CHEBFUN2 representing the maximum of F along 
%   the dimension DIM, where DIM = 1 means along the x-direction, DIM = 2 
%   for the y-direction, and DIM = 3 means along the z-direction.
%
%   WARNING: This function is not always accurate to the expected
%   precision.
% 
%   Use MAX3 for the global maximum.
%
% See also CHEBFUN3/MAX2 and CHEBFUN3/MAX3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    error('CHEBFUN:CHEBFUN3:max:input', 'CHEBFUN3 is empty');
end

% Default to max of one chebfun3:
if ( nargin < 2 )
    g = []; 
end

% Default to maximum along the x direction: 
if ( nargin < 3 )
    dim = 1;
end

% Do not allow max(F, G): 
if ( nargin > 1 && ~isempty(g) )
    error('CHEBFUN:CHEBFUN3:max:twoCHEBFUN3Inputs', ...
        'Unable to maximize two CHEBFUN3 objects.');
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
    h = chebfun2(squeeze(max(vals, [], 1)), dom(3:6));
    h = simplify(h); 
elseif ( dim == 2 )
    vals = sample(f, n, n, n);
    h = chebfun2(squeeze(max(vals, [], 2)), [dom(1:2), dom(5:6)]);
    h = simplify(h);
elseif ( dim == 3 )
    vals = sample(f, n, n, n);  
    h = chebfun2(max(vals, [], 3), dom(1:4));
    h = simplify(h);
elseif ( dim == 0 ) 
    error('CHEBFUN:CHEBFUN3:max:dim', ...
        'Dimension argument must be a positive integer scalar within indexing range.')
else
   % return the CHEBFUN3 object. This is analogous to that MAX() command in
   % MATLAB.
   h = f;  
end

end