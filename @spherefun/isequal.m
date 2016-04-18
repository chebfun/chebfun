function varargout = isequal(varargin)
%ISEQUAL Equality test for SPHEREFUN2.  
% 
% BOL = ISEQUAL(F, G) returns 0 or 1. If returns 1 then F and G are the same
% SPHEREFUN2, up to relative machine precision. If returns 0 then F and G are
% not the same up to relative machine precision. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

    [varargout{1:nargout}] = isequal@separableApprox(varargin{:});
end
