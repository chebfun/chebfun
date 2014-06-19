function spy(A, c)
%SPY    Visualize a linear CHEBOP.
%   SPY(A) creates an infographic picture of the structure of the linear CHEBOP
%   A.
%
%   SPY(A, C) is the same above, using the colour C.
%
%   To see the actual nonzero pattern of A under the current discretization, use
%   SPY(LINOP(A)).
%
% See also LINOP/SPY.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Try to create a LINOP. An error will be thrown if n is not linear.
A = linop(A);

dom = A.domain;
ends = dom;
% ends = A.domain.endsandbreaks;
de = dom(end)-dom(1);
LW = 'linewidth'; lw = 3;
C = 'color'; EC = 'EdgeColor';
if nargin < 2, c = 'b'; end

N = 3;
ish = ishold;

Ab = A.blocks;
for j = 1:size(Ab, 1)
    for k = 1:size(Ab, 2)
        isz = iszero(Ab{j, k});
        isd = ~isz && Ab{j,k}.diffOrder == 0;
        if ~A.iszero(j,k)
            if isd
                plot(ends([end 1])+(k-1)*de,-ends([end 1])-(j-1)*de, ...
                    LW, lw, C, c);
                hold on
            else
                for l = 1:numel(ends)-1
                    fill(ends([l+1 l l l+1])+(k-1)*de, ...
                        -ends([l+1 l+1 l l])-(j-1)*de, c, EC, c); 
                    hold on
                end
            end
        end
    end
end

set(gca,'xTick',[])
set(gca,'yTick',[])

if ~ish, hold off, axis equal, axis tight, end

end