function G = quasimatrix(F, varargin)
%QUASIMATRIX   A quasimatrix is an array of CHEBFUN objects.
%   F = QUASIMATRIX(OP, ...) constructs a quasimatrix representation F of
%   CHEBFUN(OP, ...).
%   
%   G = QUASIMATRIX({F1, F2, F3, ...}) constructs a quasimatrix G from the
%   cell-array of CHEBFUN objects F1, F2, ...
%
%   In Chebfun we use the term 'quasimatrix' to describe an array of CHEBFUN
%   objects. In V5 the situation becomes slightly confusing due to the
%   introduction of 'array-valued' CHEBFUN objects. To clarify:
%      quasimatrix = array of CHEBFUNs   vs   array-valued CHEBFUN
%
%   Q: Why do both of these things exist? 
%   A: For performance.
%      In many cases we would like to manipulate a collection of CHEBFUN-like
%      objects which all have the same break points and discretization sizes.
%      This is precisely what an array-valued CHEBFUN does, and at the 'tech'
%      level everything can be done using fast dense linear algebra. CHEBFUN2,
%      for example, use this type of structure in its CDR decomposition.
%
%   Q: So why bother with arrays of CHEBFUN objects (i.e., quasimatrices)?
%   A: It's not always convenient to store things in the array-valued format,
%      and sometimes it is not even possible. Currently the classes that handle
%      singularities and delta functions are not vectorized, and the support for
%      multiple columns must come at the CHEBFUN level. (This situation is not
%      ideal and we hope to rectify it in a future release.)
%
%   Q: What's the difference in how the two are used?
%   A: We've tried very hard to make user interaction with quasimatrices and
%      array-valued CHEBFUNs be as similar as possible. If you run the examples
%      below, you'll see, for instance, that both produce identical output
%      when displayed.  The easiest way to check if a CHEBFUN F, is a
%      quasimatrix is via ISQUASI(F). The size SIZE(F) of a quasimatrix and its
%      equivalent array-valued CHEBFUN will always be the same, but a
%      quasimatrix will always have NUMEL(F) > 1. You can check the number of
%      columns (or rows, if transposed) for both a quasimatrix and an
%      array-valued CHEBFUN with NUMCOLUMNS(F).
%
%   Q: Why do some methods support array-valued CHEBFUNS but not quasimatrices?
%   A: Partly because there are some operations that only make sense for the
%      array-valued representation, and partly because we haven't gotten around
%      to doing all the quasimatrix stuff yet. (Sorry!)
%
%   Q: How do I convert between the two representations?
%   A: You can use the CHEB2QUASI(F) or QUASI2CHEB(F) methods.
%
% Examples:
%   % CHEBFUNs created via the constructor are _always_ array-valued:
%   f = chebfun(@(x) [sin(x) abs(x)], 'splitting', 'on');
%   isQuasi(f) % ans = 0
%
%   % If you want to construct a quasimatrix directly, you can use this method:
%   f = quasimatrix(@(x) [sin(x) abs(x)], 'splitting', 'on');
%   isQuasi(f) % ans = 1
%
%   % HORZCAT of CHEBFUNs with the same domain gives an array-valued result:
%   x = chebfun('x');
%   f = [sin(x) cos(x)];
%   isQuasi(f) % ans = 0
%
%   % HORZCAT of CHEBFUNs with different breakpoints gives a quasimatrix:
%   x = chebfun('x');
%   f = [sin(x) abs(x)];
%   isQuasi(f) % ans = 1
%
% See also CHEBFUN, ISQUASI, CHEB2QUASI, QUASI2CHEB, NUMCOLUMNS

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER NOTE:
%   This method is basically  a wrapper for CHEB2QUASI() and CHEBFUN().
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ( iscell(F) && any(cellfun(@(f) isa(f, 'chebfun'), F)) )
    % Construct from a cell of CHEBFUNs.
    G = chebfun.cell2quasi(F);
    return
end

if ( ~isa(F, 'chebfun') )
    % Create a CHEBFUN:
    F = chebfun(F, varargin{:});
end

% Convert to a quasimatrix:
G = cheb2quasi(F);

end
