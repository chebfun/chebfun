function projOrder = getProjOrder(varargin)
%GETPROJORDER   Get projection order of a CHEBMATRIX.
%   GETPROJORDER() returns zero, ragardless of the input. It is required because
%   sublasses of CHEBDISCRETIZATION (in particular LINOP) require different
%   behaviour.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% In general, CHEBDISCRETIZATION objects are not projected.
projOrder = 0;
     
end