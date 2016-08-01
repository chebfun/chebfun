function vals = feval(varargin)
%FEVAL pointwise evaluate a DISKFUNV.
%   feval(F, t,r, 'polar') returns the evaluation of F at polar coordinates
%   (t = angular, r = radial).
%   feval(F,X, Y) returns the evaluation of F at the Cartesian
%   coordinates (X, Y)
%   See also SUBSREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( varargin{1}) ) 
    vals = []; 
    return
end

F = varargin{1}; 
x = varargin{2};
nF = F.nComponents; 
vals = zeros(nF, length(x)); 

 vals(1, :) = feval(F.components{1}, varargin{2:end});    
 vals(2, :) = feval(F.components{2}, varargin{2:end});   
end
