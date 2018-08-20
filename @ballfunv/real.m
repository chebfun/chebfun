function w = real(v)
% REAL Real part of a BALLFUNV
%   REAL(v) is the real part of the BALLFUNV v
V = v.comp;
w = ballfunv(real(V{1}),real(V{2}),real(V{3}));
end
