function N1 = plus(N1, N2) 
%PLUS plus for CHEBOP2
%
% N = PLUS(N1, N2) is the same as N = N1 + N2. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.
 
% check domains are the same. 
if all(N1.domain ~= N1.domain)
   error('CHEBOP2:PLUS','Domains of operators must be identical.'); 
end

% Plus the coefficient matrices together. 
A = N1.coeffs; B = N2.coeffs;

if iscell(A) && iscell(B) 
    C = cell(max(size(A,1),size(B,1)),max(size(A,2),size(B,2)));
    if size(B,1) < size(C,1) 
        B{size(C,1),1} = [];
    end
    if size(B,2) < size(C,2) 
        B{1,size(C,2)} = [];
    end 
    if size(A,1) < size(C,1) 
        A{size(C,1),1} = [];
    end 
    if size(A,2) < size(C,2) 
        A{1,size(C,2)} = [];
    end 
    for jj = 1:size(C,1)
        for kk = 1:size(C,2)
            if isempty(B{jj,kk})
                C{jj,kk} = A{jj,kk};
            elseif isempty(A{jj,kk})
                C{jj,kk} = B{jj,kk};
            else
                C{jj,kk} = A{jj,kk} + B{jj,kk};
            end
        end
    end
elseif iscell(A)
    BB = cell(size(B,2),size(B,1));
    for jj = 1:size(B,1)
        for kk = 1:size(B,2)
            BB(kk,jj) = {B(jj,kk)};
        end
    end
    N2.coeffs = BB;
    N1 = plus(N1,N2); 
    return;
elseif iscell(B)
    AA = cell(size(A,2),size(A,1));
    for jj = 1:size(A,1)
        for kk = 1:size(A,2)
            AA(kk,jj) = {A(jj,kk)};
        end
    end
    N1.coeffs = AA;
    N1 = plus(N1,N2); 
    return;
else
    C = zeros(max(size(A,1),size(B,1)),max(size(A,2),size(B,2)));
    if size(B,1) < size(C,1) 
        B(size(B,1)+1:size(C,1),:) = 0;
    end
    if size(B,2) < size(C,2) 
        B(:,size(B,2)+1:size(C,2)) = 0;
    end
    if size(A,1) < size(C,1) 
        A(size(A,1)+1:size(C,1),:) = 0;
    end 
    if size(A,2) < size(C,2) 
        A(:,size(A,2)+1:size(C,2)) = 0;
    end 
    C = A + B; 
end

if ~isempty(N1.lbc) || ~isempty(N1.rbc) || ~isempty(N1.ubc) || ~isempty(N1.dbc)
    warning('BCs were ignored.');
    N1.lbc =[]; N1.lbcshow = [];
    N1.rbc =[]; N1.rbcshow = [];
    N1.ubc =[]; N1.ubcshow = [];
    N1.dbc =[]; N1.dbcshow = [];
end

if ~isempty(N2.lbc) || ~isempty(N2.rbc) || ~isempty(N2.ubc) || ~isempty(N2.dbc)
    warning('BCs were ignored.');
    N2.lbc =[]; N2.lbcshow = [];
    N2.rbc =[]; N2.rbcshow = [];
    N2.ubc =[]; N2.ubcshow = [];
    N2.dbc =[]; N2.dbcshow = [];
end


N1.coeffs = C;
N1.xorder = max(N1.xorder, N2.xorder); 
N1.yorder = max(N1.yorder, N2.yorder); 

% Plus the ops together
op1 = N1.op; op2 = N2.op; 
N1.op = @(u) op1(u) + op2(u); 

end