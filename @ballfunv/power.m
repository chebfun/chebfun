function w = power(v,n)
% POWER Return a BALLFUNV to the power n
%   POWER(v) is the BALLFUNV (Vx^n, Vy^n, Vz^n)
V = v.comp;
w = ballfunv(power(V{1},n),power(V{2},n),power(V{3},n));
end
