function vals = feval(varargin)
%FEVAL  Pointwise evaluate a DISKFUNV.
%   FEVAL(F, X, Y) returns the evaluation of the DISKFUNV F at the Cartesian
%   coordinates (X, Y). 
% 
%   FEVAL(F, T, R, 'polar') returns the evaluation of F at polar coordinates
%   (T, R), where T is the angular variable and R is the radial variable. 
%
% See also SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( varargin{1}) ) 
    vals = []; 
    return
end

% Parse user inputs: 
F = varargin{1}; 
x = varargin{2};

% Extract components: 
nF = F.nComponents; 
vals = zeros(nF, length(x)); 

% Evaluate the DISKFUNV object:
vals(1, :) = feval(F.components{1}, varargin{2:end});    
vals(2, :) = feval(F.components{2}, varargin{2:end});   

end