function varargout = transpose(varargin)
% .'   DISKFUN transpose. 
%    F.' is the non-conjugate transpose of the underlying  quasimatrix
% for the diskfun F. Given F in polar coordinates, F(theta, r), 
% F.' = F(r, theta). To transpose a diskfun with respect to Cartesian
% coordinates, use the command FLIPXY. 
% 
% See also CTRANSPOSE, FLIPXY.  

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = transpose@separableApprox(varargin{:});

end
