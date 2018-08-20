function h = times(f, g)
% TIMES Component-wise multiplication of two BALLFUNV
%   TIMES(f, g) is the component-wise multiplication between the BALLFUNV 
%   f and g

%   F.*G if F is a ballfunv and G is double returns the ballfunv after
%   componentwise multiplication.F = f.comp;
F = f.comp;
G = g.comp;
h = ballfunv(F{1}*G{1}, F{2}*G{2}, F{3}*G{3});
end
