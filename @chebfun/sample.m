function varargout = sample(f, varargin)
%SAMPLE   Sample a CHEBFUN on an "appropriate" grid.
%   If F has a single FUN, VALUES = SAMPLE(F, N) returns a vector VALUES of the
%   values of F at "appropriately-chosen" points.  What "appropriately-chosen"
%   means depends on the type of representation on which F is ultimately based.
%
%   [VALUES, POINTS] = SAMPLE(F, N) returns also the vector POINTS at which the
%   values were computed.
%
%   [...] = SAMPLE(F) uses N = LENGTH(F).
%
%   If F has multiple FUNs, this function throws an error.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:  The purpose of this function (and the resulting "stack" of
% SAMPLE functions that goes down to the bottom layers) is to allow one to get
% samples of a CHEBFUN on the "right" grid in a tech-agnostic manner, i.e.,
% SAMPLE(F) should return samples on a second-kind Chebyshev grid if F is based
% on CHEBTECH2 and samples on a first-kind Chebyshev grid if F is based on
% CHEBTECH1.  The point is that you do not need to tell the function what the
% underlying tech is.
%
% The main client of this function is CHEBFUN2, which needs samples from the
% "correct" grid in a couple of places but also needs to remain tech-agnostic.
%
% This is an internal function and is not intended to be called by end-users.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ( numel(f.funs) > 1 )
        error('CHEBFUN:CHEBFUN:sample:multipleFuns', ...
              'SAMPLE is not supported for CHEBFUN objects with several FUNs.');
    end

    [varargout{1:nargout}] = sample(f.funs{1}, varargin{:});
end
