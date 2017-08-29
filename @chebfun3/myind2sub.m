function [row, col] = myind2sub(siz, ndx)
%MYIND2SUB   Townsend's version of ind2sub for linear indexing.
%
%   Notes copied from CHEBFUN2: in2sub is slow because it has a varargout. 
%   Since this is at the very inner part of the constructor and slowing 
%   things down we will make our own. This version is about 1000 times 
%   faster than MATLAB ind2sub.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

vi = rem(ndx - 1, siz(1)) + 1 ;
col = (ndx - vi) / siz(1) + 1;
row = (vi - 1) + 1;

end