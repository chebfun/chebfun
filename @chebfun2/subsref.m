function varargout = subsref(varargin)
%SUBSREF   CHEBFUN2 subsref.
% ( )
%   F(X, Y) returns the values of the CHEBFUN2 F evaluated at (X,Y). See
%   CHEBFUN/FEVAL for further details. F(:, Y) returns a chebfun 
%   representing the function F along that column slice, and F(X, :) 
%   returns a chebfun representing F along that row slice. F(:, :) returns 
%   F.
%
%   F(G) computes the composition of F and G where G is a CHEBFUN with two 
%   columns, a CHEBFUN2V, a CHEBFUN3V with two components, or a DISKFUNV. 
%   If G is a CHEBFUN with one column, a CHEBFUN2, a CHEBFUN3, a DISKFUN or
%   a SPHEREFUN, F(G) is interpreted as F(real(G), imag(G)), regardless of
%   whether G is real or complex.
%
%   F(X, Y) with CHEBFUNs X and Y returns the CHEBFUN G(t) = F(X(t), Y(t)).
%   If X and Y are CHEBFUN2 objects, then F(X, Y) is a CHEBFUN2.
%   If X and Y are CHEBFUN3 objects, then F(X, Y) is a CHEBFUN3.
%
% .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
% {}
%   F{S1, S2, S3, S4} restricts F to the domain [S1, S2, S3, S4]. See
%   CHEBFUN2/RESTRICT for further details. Note that F{[S1,S2, S3, S4]} is 
%   not supported due to the behaviour of the MATLAB subsref() command.
%
% See also FEVAL, GET, RESTRICT, SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = subsref@separableApprox(varargin{:});

end