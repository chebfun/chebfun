function varargout = subsref(varargin)
%SUBSREF       CHEBFUN2 subsref.
% ( )
%   F(X, Y) returns the values of the CHEBFUN2 F evaluated at (X,Y). See
%   CHEBFUN/FEVAL for further details. F(:, Y) returns a chebfun representing
%   the function F along that column slice, and F(X, :) returns a chebfun
%   representing F along that row slice. F(:, :) returns F.
%
%   F(G), where G is also a CHEBFUN2V with two components
%   computes the composition of F and G.
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2, S3, S4} restricts F to the domain [S1, S2, S3, S4]. See
%   CHEBFUN2/RESTRICT for further details. Note that F{[S1,S2, S3, S4]} is not
%   supported due to the behaviour of the MATLAB subsref() command.
%
% See also FEVAL, GET, RESTRICT, SUBSREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = subsref@separableApprox(varargin{:});

end
