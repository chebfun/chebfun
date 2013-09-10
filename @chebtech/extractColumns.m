function f = extractColumns(f, colIdx)

% TODO: Document

f.values = f.values(:, colIdx);
f.coeffs = f.coeffs(:, colIdx);
f.vscale = f.vscale(colIdx);

end