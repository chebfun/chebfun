function F = diff(F, dim, n)
%DIFF   Componentwise derivative of a DISKFUNV.
%   DIFF(F) is the derivative of each component of F in
%   x-direction.
%
%   DIFF(F, DIM) is the first derivative of F along the
%   dimension DIM.
%     DIM = 1 (default) is the derivative in the x-direction.
%     DIM = 2 is the derivative in the y-direction.
%    
%
%   DIFF(F, DIM, N) is the Nth derivative each component of F
%   along the dimension specified.
%
%   See also CURL and DIV

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check:
if ( isempty(F) ) 
    return
end

% Defaults:
if ( nargin == 1 || isempty(dim) )
    dim = 1;
end
if ( nargin < 3 ) 
    n = 1; 
end

% Diff each component. 
F.components{1} = diff(F.components{1}, dim, n); 
F.components{2} = diff(F.components{2}, dim, n); 

end
