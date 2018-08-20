function pass = test_normalized_Legendre()
% Test the equality for P^m_l, 0 <= l <= n and -m <= l <= m
eps = 1e-10;

n = 10;
A = ballfun.normalized_legendre(n);

Max_difference = 0;

for l = 0:n
    for m = 0:l
        Z = convert_normalized_Legendre(m,A(:,l*(l+1)/2+m+1));
        Y = spherefun.sphharm(l,m);
        % Divide by sqrt(2) if m > 0
        if m > 0
            Y = Y/sqrt(2);
        end
        Max_difference = max(norm(Z-Y),Max_difference);
        %fprintf("l = %d, m = %d : norm = %g\n", l,m,norm(Y-Z))
    end
end
pass = Max_difference < eps;
end

function Z = convert_normalized_Legendre(m,f)
p_tilde = size(f,1);
F = zeros(p_tilde,2*m+1);
F(:,2*m+1) = f;
Z = spherefun.coeffs2spherefun(F);
end
