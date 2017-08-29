function varargout = resampling(varargin) %#ok<STOUT>
%RESAMPLING   CHEBFUN 'resampling' option.
%   RESAMPLING ON and RESMAPLING OFF are no longer supported. See the section on
%   refinementFunction in CHEBTECH.TECHPREF documentation for further details.
%
% See also CHEBTECH.TECHPREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:resampling:deprecated', ...
    ['The syntax ''resampling on'' is no longer supported.\n', ...
     'See the section on refinementFunction in ''help chebtech.techPref'' for details.']);
 
end


