function out = vertcat(varargin)
%VERTCAT   Vertical concatenation of CHEBFUN objects.
%   VERTCAT(F1, F2, ...) of column CHEBFUNs F1, F2, ... produces a CHEBMATRIX.
%
%   VERTCAT(F1, F2, ...) of row CHEBFUNs F1, F2, ... produces a quasimatrix.
%
% See alos HORZCAT, CAT, CHEBMATRIX, QUASIMATRIX.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%  Chebmatrices may only have 'single' (i.e., scalar) entries in their blocks.
%  NUM2CELL should be enough to do this, but it not yet extensively tested.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Test.

% Find the locations of the CHEBFUN objects in the inputs:
isCheb = cellfun('isclass', varargin, 'chebfun');
chebfun1 = varargin{find(isCheb, 1, 'first')};

% Check transpose state:
if ( ~all(chebfun1(1).isTransposed == ...
        cellfun(@(f) double(f(1).isTransposed), varargin(isCheb)) ) )
    error('CHEBFUN:CHEBFUN:vertcat:transpose', ...
        'Dimensions of matrices being concatenated are not consistent. ');
end

numElements = cellfun(@(u) min(size(u)), varargin);
if ( any(numElements > 1) )
    for k = 1:numel(varargin)
        varargin{k} = chebmatrix(num2cell(varargin{k}));
    end
    out = vertcat(varargin{:});
    % TODO: This warning should be removed eventually.
    warning('CHEBFUN:CHEBFUN:vertcat:join', ...
        ['Vertical concatenation of CHEBFUN objects now produces a CHEBMATRIX\n', ...
         'The V4 behaviour can be reproduced using the JOIN() method.']);
    warning('off', 'CHEBFUN:CHEBFUN:vertcat:join');
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