function w = imag(v)
% IMAG Imaginary part of a BALLFUNV
%   IMAG(v) is the imaginary part of the BALLFUNV v

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

V = v.comp;
w = ballfunv(imag(V{1}),imag(V{2}),imag(V{3}));
end
