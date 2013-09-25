function L = rbc(L,f,value)
d = domain(L);
E = linop.evalAt(d(end),d);
if nargin < 3
    value = f;
    f = E;
else
    f = E*f;
end
L = L.bc(f,value);
end
