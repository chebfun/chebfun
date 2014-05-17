function varargout = fracCumSum(varargin)%#ok<STOUT>
%FRACUMSUM   FRACCUMSUM does not support UNBNDFUN objects.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

error('CHEBFUN:UNBNDFUN:fracCumSum:nosupport', ...
      'FRACCUMSUM does not support UNBNDFUN objects.')

end