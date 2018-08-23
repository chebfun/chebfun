function varargout = subsref(f, index)
%SUBSREF   BALLFUN subsref.
%( )
%   F(X, Y, Z) or F(X, Y, Z, 'cart') returns the values of the BALLFUN 
%   object F evaluated at the points (X, Y, Z) in cartesian coordinates.
%
%   F(R, L, TH, 'polar') returns the values of the BALLFUN object F 
%   evaluated at the points (R, L, TH) in spherical scoordinates.
%
%   F(R, :, :) returns a spherefun representing the function F along a 
%   radial shell. 
% 
%   F(:, :, :) returns F.
%
%   F(G) where G is a BALLFUN returns the BALLFUN representing the
%   composition F(G). 
%
%{ } 
%   Not supported.
%
%   F.PROP returns the property of F specified in PROP.
%
% See also DISKFUN/FEVAL, BALLFUN/GET. 

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Implemented this. 

end