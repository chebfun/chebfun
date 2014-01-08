function out = vertcat(varargin)
%VERTCAT   Vertical concatenation of CHEBFUN objects.
%   VERTCAT of a CHEBFUN is not yet supported.

% TODO: Document.
% TODO: Test.

% Find the locations of the CHEBFUN objects in the inputs:
chebfunLocs = cellfun('isclass', varargin, 'chebfun');
chebfun1 = varargin{find(chebfunLocs, 1, 'first')};

% Horizontal concatenation of row CHEBFUN objects produces a CHEBMATRIX:
if ( chebfun1(1).isTransposed )
    args = cellfun(@transpose, varargin, 'UniformOutput', false);
    out = horzcat(args{:}).';
else
    out = chebmatrix(varargin.');
end

end
