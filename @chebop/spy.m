function spy(A)
%SPY    Visualize a linear CHEBOP.
%   SPY(A) creates an infographic picture of the structure of the linear CHEBOP
%   A.
%
%   To see the actual nonzero pattern of A under the current discretization, use
%   SPY(LINOP(A)).
%
% See also LINOP/SPY.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Try to create a LINOP. An error will be thrown if n is not linear.
A = linop(A);

% Domain information:
dom = A.domain;
domDiff = dom(end) - dom(1);

% Were we holding when we entered this method?
ish = ishold;

% Begin by looping through the DE blocks:
Ab = A.blocks;
for j = 1:size(Ab, 1)
    for k = 1:size(Ab, 2)
        % Is the {j,k} block a zero operator?
        isZero = iszero(Ab{j, k});
        % Is the {j,k} block a multiplication operator?
        isMult = (~isZero && ( Ab{j,k}.diffOrder == 0) );
        if ~isZero
            if isMult
                % Draw a diagonal line for multiplication operators.
                plot(dom([end 1]) + (k - 1)*domDiff, ...
                    -dom([end 1]) - (j - 1)*domDiff, ...
                    'Linewidth', 3, 'Color', 'b');
                hold on
            else
                % Drow a filled rectangle for other operators.
                for l = 1:numel(dom)-1
                    fill(dom([l+1 l l l+1]) + (k - 1)*domDiff, ...
                        -dom([l+1 l+1 l l]) - (j - 1)*domDiff, ...
                        'b', 'EdgeColor', 'b'); 
                    hold on
                end
            end
        end
    end
end

% TODO: Plot lines for constraints?

% Remove the ticks from the axes
set(gca,'xTick',[])
set(gca,'yTick',[])

% If we were not holding, make the plot look nice
if ~ish
    hold off
    axis equal
    axis tight
end

end
