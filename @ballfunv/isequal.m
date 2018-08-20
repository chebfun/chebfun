function b=isequal(f, g)
% ISEQUAL Test the equality between two BALLFUNV
%   ISEQUAL(f, g) is the boolean f == g
b = iszero(f-g);
end
