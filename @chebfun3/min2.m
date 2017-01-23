function h = min2(f, g, dims)
%MIN2   Minimum value of a CHEBFUN3 in two directions.
%   MIN2(F) returns a 1D chebfun representing the minimum of the CHEBFUN3 
%   along the x and y directions, i.e, MIN2(F) = @(z) min(F(:, :, z)).
%
%   MIN2(F, [], DIMS) returns a CHEBFUN representing the minimum of F along 
%   the DIMS directions. DIMS = [1, 2] means along the x and y directions, 
%   etc.
%
%   WARNING: This function is not always accurate to the expected precision. 
% 
%   Use MIN3 for computing the global minimum.
%
% See also CHEBFUN3/MIN and CHEBFUN3/MIN3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    error('CHEBFUN:CHEBFUN3:min2:input', 'CHEBFUN3 is empty');
end

% Default to max2 of one chebfun3:
if ( nargin < 2 )
    g = []; 
end

% Default to maximum along the x and y directions:
if ( nargin < 3 )
    dims = [1, 2];
end

% Do not allow max(F, G): 
if ( nargin > 1 && ~isempty(g) )
    error('CHEBFUN:CHEBFUN3:max2:twoCHEBFUN3Inputs', ...
        'Unable to maximize two CHEBFUN3 objects.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We have no idea how to achieve this in an efficient way. This
% is an attempt to return a result, but typically it won't be accurate to 
% more than 4-5 digits. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dom = f.domain;
n = 129;
if ( all(dims == [1, 2]) || all(dims == [2, 1]) )
    vals = sample(f, n, n, n);
    temp = chebfun3.unfold(vals, [3]);
    h = chebfun(min(temp, [], 2), dom(5:6), 'splitting', 'on');
    h = simplify(h); 
elseif ( all(dims == [1,3]) || all(dims == [3, 1]) )
    vals = sample(f, n, n, n);
    temp = chebfun3.unfold(vals, [2]);
    h = chebfun(min(temp, [], 2), dom(3:4), 'splitting', 'on');
    h = simplify(h);
elseif ( all(dims == [2, 3]) || all(dims == [3, 2]) )
    vals = sample(f, n, n, n);
    temp = chebfun3.unfold(vals, [1]);
    h = chebfun(min(temp, [], 2), dom(1:2), 'splitting', 'on');
    h = simplify(h);
elseif ( dims == 0 )
    error('CHEBFUN:CHEBFUN3:min2:dims', ...
        'Dimension arguments must be two positive integer scalars within indexing range.')
else
   % return the CHEBFUN3. This is analogous to that MIN() command in
   % MATLAB.
   h = f;  
end

end
