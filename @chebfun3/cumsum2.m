function f = cumsum2(f, dims)
%CUMSUM2   Double indefinite integral of a CHEBFUN3.
%   F = CUMSUM2(F) returns the double indefinite integral of a CHEBFUN3 F. 
%   By default, that means cumsum in the first 2 variables, i.e., x and y:
%                   y  x
%                  /  /
%   CUMSUM2(F) =  |  |   F(x,y,z) dx dy
%                 /  /
%                c  a
%
%   where [a,b] x [c,d] x [e,g] is the domain of F.
% 
%   DIMS is a vector containing two of the three indices 1,2,3 to show
%   which two of the dimensions are to be used.
% 
% See also CHEBFUN3/CUMSUM, CHEBFUN3/CUMSUM3, CHEBFUN3/SUM, CHEBFUN3/SUM2 
% and CHEBFUN3/SUM3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) ) 
    f = [];
    return
end

% Default to cumsum in the first two variables x and y:
if ( nargin == 1 )
    dims = [1, 2];
end

if ( numel(dims) ~= 2 )
    error('CHEBFUN:CHEBFUN3:cumsum2:dims', 'Dims should have 2 entries.');
end

if ismember(1, dims) 
    % cumsum along the 1st variable
    f.cols = cumsum(f.cols);
end
if ismember(2, dims)
    % cumsum along the 2nd variable
    f.rows = cumsum(f.rows);
end
if ismember(3, dims) 
    % cumsum along the 3rd variable
    f.tubes = cumsum(f.tubes);
end

end