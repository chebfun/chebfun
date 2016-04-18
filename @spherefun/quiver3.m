function varargout = quiver3(varargin)
%QUIVER3  3-D quiver plot of a CHEBFUN2V at data mapped by a SPHEREFUN2.
%
% QUIVER3(Z, F) plots velocity vectors at the equally spaced surface point
% specified by the SPHEREFUN2 Z. We use Z to map a uniform grid. F should be
% a CHEBFUN2V.
%
% QUIVER3(X, Y, Z, F) plots velocity vectors at (x,y,z), where X, Y, Z are
% SPHEREFUN2 objects which we use to to map a uniform grid. F should be a
% CHEBFUN2V.
%
% Alternative syntax for this command is:
% QUIVER3(X,Y,Z,[f;g;h]) or QUIVER3(X,Y,Z,f,g,h), where f, g, and h are
% SPHEREFUN2 objects.
%
% QUIVER(...,'numpts',N) plots arrows on a N x N uniform grid.
%
% This command is a wrapper to CHEBFUN2V/QUIVER3, and is required because
% SPHEREFUN2 methods take priority over CHEBFUN2V methods.
%
% See also CHEBFUN2V/QUIVER3.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = quiver3@separableApprox(varargin{:});
end
