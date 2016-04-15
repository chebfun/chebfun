function [row, col] = myind2sub(siz, ndx)
% Copied from Chebfun2, written by Alex Townsend.
% My version of ind2sub. In2sub is slow because it has a varargout. Since this
% is at the very inner part of the constructor and slowing things down we will
% make our own. This version is about 1000 times faster than MATLAB ind2sub.

vi = rem(ndx - 1, siz(1)) + 1 ;
col = (ndx - vi) / siz(1) + 1;
row = (vi - 1) + 1;

end