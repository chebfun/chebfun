function vals = feval(varargin)
%FEVAL pointwise evaluate a DISKFUN.
%   feval(F, X,Y) returns the evaluation of F at the polar coordinates (theta,r).
%   F(F,X, Y, 'cart') returns the evaluation of F at the Cartesian
%   coordinates (X, Y)
% See also SUBSREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( varargin{1}) ) 
    vals = []; 
    return
end

%figure out if cartesian or polar
 %iscart = diskfun.coordsetting(varargin{:});

F = varargin{1}; 
x = varargin{2};
y = varargin{3};
nF = F.nComponents; 
vals = zeros(nF, length(x)); 
coordsetting = varargin{4};
% Evaluate each component:

%if nargin < 4
 %   coords = 1;
%end
%if strcmpi(coords,'cart')
%    coords = 0; 
%end
%if iscart                  
 %   coords = 'cart';
%else
 %   coords = 'polar';
%end
    for jj = 1:nF
     vals(jj, :) = feval(F.components{jj}, x, y, varargin{3:end});  
    end

end
