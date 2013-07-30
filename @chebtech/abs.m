function f = abs(f)

f.values = abs(f.values);
f.coeffs = f.chebpoly(f.values);

end