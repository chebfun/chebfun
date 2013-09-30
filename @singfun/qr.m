function [f, R, E] = qr(f, outputFlag, methodFlag)
%QR   SINGFUN does not support array-valued objects.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('CHEBFUN:SINGFUN:qr:nosupport', ...
      'QR does not support SINGFUN objects.')
    
end