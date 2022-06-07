function ct2 = createCT2(W,data)
data.vscale = max(abs(W(:)));
ct2 = chebtech2(W, data);
ct2.coeffs = sum(abs(ct2.coeffs), 2);
end