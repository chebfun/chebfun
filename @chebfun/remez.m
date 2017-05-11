function varargout = remez(f, varargin)
%REMEZ   Best polynomial or rational approximation for real valued CHEBFUNs.
%   This code is obsolete as of May 2017.  Use MINIMAX instead.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

warning('CHEBFUN:CHEBFUN:remez', ...
        [' This command is deprecated.', ...
         ' Use minimax instead.']);
polyOutput = detectType(f,varargin{:});

if ( polyOutput )
    [p,err,status] = minimax(f,varargin{:});
    varargout = {p, err, status};
else
    [p,q,r,err,status] = minimax(f,varargin{:});
    varargout = {p, q, r, err, status};
end

end
    
function polyOutput = detectType(~, varargin)

isSilent = 0;
for k = 1:length(varargin)
    if ( ischar(varargin{k}) && strcmpi('silent', varargin{k}) )
        isSilent = 1;
    end
end

% Detect polynomial / rational approximation type.
polyOutput = true;
if ( mod(nargin - isSilent, 2) ) % Odd number of inputs --> rational case.                             
    polyOutput = false;
end

end
