function M = mult(A, f, lambda)

n = sum(dim(A));

% get Chebyshev T coefficients
a = flipud(get(f, 'coeffs'));

if ( numel(a) == 1 )
    M = a*speye(n);
    return
end

% prolong or truncate coefficients
if ( numel(a) < n )
    a = [a ; zeros(n - numel(a), 1)];
else
    a = a(1:n);  % truncate.
end

if ( lambda == 0 )
    a = a/2;  % just to make formula easier.
    M = sptoeplitz([2*a(1);a(2:end)], [2*a(1);a(2:end)]);
    H = sphankel(a(2:end));
    sub1 = 2:length(a); sub2 = 1:length(a)-1;
    M(sub1, sub2) = M(sub1, sub2)+ H;
elseif ( lambda == 1 )
    M = sptoeplitz([2*a(1);a(2:end)], [2*a(1);a(2:end)])/2;
    sub = 1:length(a)-2;
    M(sub, sub) = M(sub, sub) - sphankel(a(3:end)/2);
else
    % Want the C^{lam}C^{lam} Cheb Multiplication matrix.
    
    dummy = blockUS(n, [-1 1]);
    
    % Convert ChebT of a to ChebC^{lam}
    a = convert(dummy, 1, lambda+1) * a;
    
    nnza = find(abs(a)>eps); band = nnza(end);
    M = 0*speye(n);
    for d = 0:n-1
        for k=max(0,d-band):min(n-1,d+band)
            if(k<=d)
                as=0;
                c = ChebC(lambda,k,d-k,0);
                for s=0:k
                    if(2*s+d-k<n)
                        as = as + a(2*s+d-k+1)*c;
                        c = c*factorcheb(lambda,k,2*s+d-k,s);
                    end
                end
                M(d+1,k+1) =  as;
            elseif(k>d)
                as=0;
                c = ChebC(lambda,k,k-d,k-d);
                for s=k-d:k
                    if(2*s+d-k<n)
                        as = as + a(2*s+d-k+1)*c;
                        c = c*factorcheb(lambda,k,2*s+d-k,s);
                    end
                end
                M(d+1,k+1) = as;
            end
        end
    end
    M = M(:,1:n);
end
end

function c3 = ChebC(v,m,n,s)
%algebraic mess:
c3 = ((m+n+v-2*s)/(m+n+v-s))*prod(((v:v+s-1)./(1:s)).*((2*v+m+n-2*s:2*v+m+n-s-1)./(v+m+n-2*s:v+m+n-s-1)));
c3 = c3*prod(((v:v+m-s-1)./(1:m-s)).*((n-s+1:m+n-2*s)./(v+n-s:v+m+n-2*s-1)));
end

function fac = factorcheb(v,i,j,s)
fac = ((v+s)/(s+1))*((i-s)/(v+i-s-1))*((2*v+i+j-s)/(v+i+j-s));
fac = fac*((v+j-s)/(j-s+1))*((i+j+v-s)/(i+j+v-s+1));

end