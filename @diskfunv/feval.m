function vals = feval(F, x, y, coords)
%FEVAL pointwise evaluate a DISKFUN.
%   feval(F, X,Y) returns the evaluation of F at the polar coordinates (theta,r).
%   F(F,X, Y, 'cart') returns the evaluation of F at the Cartesian
%   coordinates (X, Y)
% See also SUBSREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) ) 
    vals = []; 
    return
end

nF = F.nComponents; 
vals = zeros(nF, length(x)); 

% Evaluate each component:

if nargin < 4
    coords = 1;
end
if strcmpi(coords,'cart')
    coords = 0; 
end

    for jj = 1:nF
     vals(jj, :) = feval(F.components{jj}, x, y, coords);  
    end


end
