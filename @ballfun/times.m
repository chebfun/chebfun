function h=times(f,g)
% TIMES Multiplication of two BALLFUN functions
%   TIMES(f, g) is the multiplication between the BALLFUN functions f
%   and g
if (nnz(size(f)-size(g))==0)
    F = f.coeffs;
    G = g.coeffs;
    Fvals = ballfun.coeffs2vals(F);
    Gvals = ballfun.coeffs2vals(G);
    Hvals = Fvals.*Gvals;
    H = ballfun.vals2coeffs(Hvals);
    h = ballfun(H);
else
    error('BALLFUN:isequal:unknown', ...
    ['Undefined function ''times'' for different size of ballfun functions : ' ...
     '%s and %s.'], mat2str(size(f)), mat2str(size(g)));
end

end
