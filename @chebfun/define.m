function f = define(f, s ,v) %#ok<INUSD>
%DEFINE   Deprecated function.
%   DEFINE(...) is deprecated and has been replaced by DEFINEPOINT() and
%   DEFINEINTERVAL(). However, most users will interact with both of these
%   methods via CHEBFUN/SUBSREF.
%
% See also DEFINEPOINT, DEFINEINTERVAL, CHEBFUN/SUBSREF.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:CHEBFUN:define:deprecated', ...
    'CHEBFUN/DEFINE is deprecated. Use DEFINEPOINT or DEFINEINTERVAL.');

end
