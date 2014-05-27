function out = vertcat(varargin)
%VERTCAT   Vertical concatenation of CHEBFUN objects.
%   VERTCAT of a CHEBFUN is not yet supported.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Document.
% TODO: Test.

% Find the locations of the CHEBFUN objects in the inputs:
chebfunLocs = cellfun('isclass', varargin, 'chebfun');
chebfun1 = varargin{find(chebfunLocs, 1, 'first')};

numElements = cellfun(@(u) min(size(u)), varargin);
if ( any(numElements > 1) )
    args = {};
    for k = 1:numel(varargin)
        varargin{k} = chebmatrix(num2cell(varargin{k}));
    end
    out = vertcat(varargin{:});
    return
end

% Horizontal concatenation of row CHEBFUN objects produces a CHEBMATRIX:
if ( chebfun1(1).isTransposed )
    args = cellfun(@transpose, varargin, 'UniformOutput', false);
    out = horzcat(args{:}).';
else
    out = chebmatrix(varargin.');
end

end