function A = convertToCorrectTech(A, newtech)
%CONVERTTOCORRECTTECH   Convert the entries of a CHEBMATRIX to another 
%TECH.
%   CONVERTTOCORRECTTECH(A, NEWTECH) converts the entries of a CHEBMATRIX A
%   to another TECH NEWTECH.
%
% See also CHEBFUN/CONVERTTOCORRECTTECH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty case.
if ( isempty(A) )
    return
end

% CONVERTTOCORRECTECH of a CHEBMATRIX with inf x inf block(s) is not 
% supported.
s = cellfun(@(b) min(size(b)), A.blocks);
if ( ~all(isfinite(s(:))) )
    error('CHEBFUN:CHEBMATRIX:convertToCorrectTech:notSupported', ...
       ['CONVERTTOCORRECTTECH of a CHEBMATRIX with inf x inf block(s)', ... 
        ' is not supported.']);
end

% If no TECH specified, do nothing.
if ( nargin == 1 )
    return
end

% Convert if necessary, using CHEBFUN/CONVERTTOCORRECTTECH.
A.blocks = cellfun(@(v) chebfun.convertToCorrectTech(v, newtech), ...
	A.blocks, 'uniformOutput', false);

end
