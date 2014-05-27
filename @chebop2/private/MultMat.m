function M = MultMat(a,bn,varargin)
% MULTMAT(A,BN) This forms a truncation of the Multiplication operator for
% a(x) in Cheb T series.  It truncates so that it works for vectors of length BN. It passes
% back a rectangular matrix. 
% MULTMAT(A,BN,lam) This forms a truncation of the Multiplication operator
% for a(x) in C^{lam} series. 

% not formed. 
    % Computes the multiplication matrix, staying in Chebyshev T. Compute
    % coefficients for a(x)*ChebTfunc
    
    if max(size(a)) == 1 
        % scalar multiplication 
        M = a; 
        return; 
    end
    
    
    
    a = a(:); 
    an = length(a);
    %remove trailing zeros
    atrun = truncate(a,eps);
    at = length(atrun); 
    
    if at == 1 
       M = atrun; 
       return;
    end
    
    %prolong coeffs. 
    %a(length(a)+1:(an+bn)) = 0;
    if(an<bn)
        a(length(a)+1:bn)=0;
    end
    if(an>bn)
        a = a(1:bn);
    end
    
    
    a = a(:);
    if(nargin>2), lam = varargin{1}; end
    
    if(nargin==2||lam==0)
    
%         if(bn<=1e2)   
%         %compute multiplication. 
%         if(length(a)>1)
%             a=.5*a;
%             r=[2*a(1);a(2:end)];
%             M = mytoeplitz(r,r);
%             a = truncate(a,eps);
%             H = hankel(a(2:end));
%             %M(2:end,1:end-1) = M(2:end,1:end-1) + H;
%             M(2:length(a),1:length(a)-1) = M(2:length(a),1:length(a)-1)+ H;
%         else
%             M = toeplitz(2*a);
%         end
%         %half because of J map. 
%         %M = M(:,1:bn);
%         %MultMat = MultMat(:,1:bn);
%         else
    
        %compute multiplication. 
        if(at>1),
            %a=.5*atrun;
            a = .5*a;
            M = sptoeplitz([2*a(1);a(2:end)],[2*a(1);a(2:end)]);
            %Construct hankel by hand for speed. 
            H = sphankel(a(2:end));
            M(2:length(a),1:length(a)-1) = M(2:length(a),1:length(a)-1)+ H;
        else
            M = sptoeplitz(2*atrun);
        end
        %half because of J map. 
        %M = .5*M(:,1:bn);
%         end
    elseif(lam==1)
        % Want the U*U Cheb Multiplication matrix.
        % Convert ChebT of a to ChebU 
        a = transMat(length(a),0)*a;
        er = flipud(cumsum(flipud(a(1:2:end))));
        or = flipud(cumsum(flipud(a(2:2:end))));

        r(1:2:length(a)) = er; r(2:2:length(a)) = or; 

        c = [r(3:end) 0 0];
        M = sptoeplitz(r,r) - sphankel(c); 
        M = M(:,1:bn);
    elseif(lam>1)
        % Want the C^{lam}C^{lam} Cheb Multiplication matrix.
        
        % Convert ChebT of a to ChebC^{lam}
        for ind=0:lam-1
            a = transMat(length(a),ind)*a;
        end

        nnza = find(abs(a)>1e-16);
        band = nnza(end);
        M = spalloc(length(a),bn,2*band*bn);
        n=length(a);
        for d=0:n-1
            for k=max(0,d-band):min(n-1,d+band)
                if(k<=d)
                    as=0;
                    c = ChebC(lam,k,d-k,0);
                    for s=0:k
                        if(2*s+d-k<n)
                            as = as + a(2*s+d-k+1)*c;
                            c = c*factorcheb(lam,k,2*s+d-k,s);
                        end
                    end
                    M(d+1,k+1) =  as;
                elseif(k>d)
                   as=0;
                   c = ChebC(lam,k,k-d,k-d);
                   for s=k-d:k
                       if(2*s+d-k<n)
                            as = as + a(2*s+d-k+1)*c;
                            c = c*factorcheb(lam,k,2*s+d-k,s);
                       end
                   end
                   M(d+1,k+1) = as;
                end
            end
        end
        M = M(:,1:bn);
    end 
end




function T = transMat(n,lam)
% TRANSMAT(N,LAM) This computes the truncation of the operator that 
% transforms C^{lam} (Ultraspherical polynomials) to C^{lam+1}.  The 
% truncation gives back a matrix of size n x n. 

if(n==1)
    T = 1;
    return;
end

%Relation is C_n^(lam) = (lam/(n+lam))(C_n^(lam+1) - C_{n-2}^(lam+1))

if(n<1e2)
    if(lam~=0)
        dg = lam./(lam + (2:n-1))';
        T = diag([1;lam./(lam+1);dg],0) + diag(-dg,2);
    elseif(lam==0)
        dg = .5*ones(n-2,1);
        T = diag([1;.5;dg],0) + diag(-dg,2);
    end
else
if(lam ~= 0) 
    dg = lam./(lam + (2:n-1))'; 
    B(:,1) = [1;lam./(lam+1);dg]; B(:,2)=[0;0;-dg];
    T = spdiags(B,[0,2],n,n);
elseif(lam==0)
    % Cheb T is special case because of different scaling.
    dg = .5*ones(n-2,1);
    B(:,1) = [1;.5;dg]; B(:,2)=[0;0;-dg];
    T = spdiags(B,[0,2],n,n);
end
end
end 





function c3 = ChebC(v,m,n,s)
        
        %algebraic mess:
        c3 = ((m+n+v-2*s)/(m+n+v-s))*prod(((v:v+s-1)./(1:s)).*((2*v+m+n-2*s:2*v+m+n-s-1)./(v+m+n-2*s:v+m+n-s-1)));
        c3 = c3*prod(((v:v+m-s-1)./(1:m-s)).*((n-s+1:m+n-2*s)./(v+n-s:v+m+n-2*s-1)));       
end

function fac = factorcheb(v,i,j,s)
%     %c_{s+1}(i,j) = fac*c_{s}(i,j), hopefully can be used to speed up things.
%     fac = ((v+s)/(s+1))*((i-s)/(v+i-s-1))*((2*v+i+j-2*s-2)/(v+i+j-2*s-2));
%     fac = fac*((2*v+i+j-2*s-1)./(v+i+j-2*s-1))*((v+i+j-s-1)/(2*v+i+j-s-1));
%     fac = fac*((j-s)/(v+j-s-1))*((v+i+j-2*s-2)/(i+j-2*s-1))*((v+i+j-2*s-1)/(i+j-2*s));
%     fac = fac*((i+j+v-2*s-2)/(i+j+v-s-1))*((i+j+v-s)/(i+j+v-2*s));
    
      %c_{s+1}(i,j) = fac*c_{s}(i,j), hopefully can be used to speed up things.
      fac = ((v+s)/(s+1))*((i-s)/(v+i-s-1))*((2*v+i+j-s)/(v+i+j-s));
      fac = fac*((v+j-s)/(j-s+1))*((i+j+v-s)/(i+j+v-s+1));
    
end


function T = sptoeplitz(col,row)
% SPTOEPLITZ Sparse Toeplitz matrix.
%    SPTOEPLITZ(C,R) produces a sparse nonsymmetric Toeplitz matrix having
%    C as its first column and R as its first row. Neither C nor R needs to
%    be sparse. No full-size dense matrices are formed.
%
%    SPTOEPLITZ(R) is a sparse symmetric/Hermitian Toeplitz matrix.
%
%    Examples:
%       sptoeplitz( real( (1i).^(0:8) ) )   % 9x9, 41 nonzeros
%       sptoeplitz( [-2 1 zeros(1,9998)] ); % classic 2nd difference
%
%    See also TOEPLITZ, SPDIAGS.

%    Copyright (c) 2006 by Tobin Driscoll (tobin.driscoll@gmail.com).
%    First version, 11 December 2006.

% This part is borrowed from built-in Toeplitz.
if nargin < 2  % symmetric case
  col(1) = conj(col(1)); row = col; col = conj(col); 
else
  if col(1)~=row(1)
    warning('MATLAB:sptoeplitz:DiagonalConflict',['First element of ' ...
      'input column does not match first element of input row. ' ...
      '\n         Column wins diagonal conflict.'])
  end
end


% Size of result.
m = length(col(:));  n = length(row(:));
%m=n;
%Only use toeplitz if you have too... fairly slow. 

if((m<2e3 && n<2e3) )%&& length(find(col))>m/2 && length(find(row))>n/2)
    if(nnz(col)==1)
        Ic=find(col);
        if(Ic==1), T = spdiags(col(Ic)*ones(m,1),0,m,n);
        else
            T = spdiags([col(Ic)*ones(m,1) col(Ic)*ones(m,1)],[-Ic+1,Ic-1],m,n);
        end
    else
        T = toeplitz(col,row);
        T = sparse(T); 
    end
else
    % Locate the nonzero diagonals.
    [ic,jc,sc] = find(col(:));
    row(1) = 0;  % not used
    [ir,jr,sr] = find(row(:));

    % Use spdiags for construction.
    d = [ ir-1; 1-ic ];
    B = repmat( [ sr; sc ].', min(m,n),1 );
    T = spdiags( B,d,m,n );
end
end



function H = sphankel(r)
% SPHANKEL(R) this forms a sparse hankel matrix by forming it as an upside-
% down toeplitz matrix. 

%Hankel is an upside-down toeplitz matrix. 
r = flipud(r(:));   %ensure column vector. 
H = fliplr(triu(sptoeplitz(r,r)));
end