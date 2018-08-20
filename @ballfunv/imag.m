function w = imag(v)
% IMAG Imaginary part of a BALLFUNV
%   IMAG(v) is the imaginary part of the BALLFUNV v
V = v.comp;
w = ballfunv(imag(V{1}),imag(V{2}),imag(V{3}));
end
