%% Interactive Interpolation
% Nick Hale, 21st November 2012

%%
% (Chebfun example approx/InteractiveInterp.m) 
% [Tags: #interpolation, #INTERP1, #LEBESGUE]

function InteractiveInterp

%%
% This example is based upon some code that Rodrigo Platte once gave me to
% demostrate what makes a good set of interpolation points (and of course,
% what doesn't!).

%%
% The game is to select interpolation points by clicking on the function.
% The code will then compute the interpolant through these points and plot
% it on the same graph. For those interested in the details, it computes
% the interpolant in stable way using barycentric form via CHEBFUN/INTERP1,
% which uses BARY and BARY_WEIGHTS

%%
% Now it's your turn! If you're reading this example, you probably already
% know what kind of points to choose, but it's still fun to play!
    
plot(0)
% ClickToInterpolate() % Commented out to ensure the Example publishes OK.

%%

    function ClickToInterpolate(F)
        % Select interpolation points by right clicking.
        % When you get bored, left or double click to finish.
        if ( nargin == 0 ) 
          F = @(x) 1-.9*abs(x);          % default function 
        end
        
        % initialise
        x = []; d = domain(-1,1);
        LW = 'LineWidth'; MS = 'MarkerSize'; FS = 'FontSize';
        % loop 
        while 1                          % keep clicking!
          hold off
          plot(xx,F(xx),'-k',LW,2), shg  % plot function F
          hold on, axis([-1 1 -1 2])
          if ( ~isempty(x) )
            plot(x,F(x),'.b',MS,20)      % plot interpolation points
            plot(x,0*x,'.k',MS,6)        % plot x values alone
            y = interp1(x,F(x),d);       % interpolate the data
            plot(y,'-b',LW,2), shg       % plot interpolant
            if ( numel(x) > 1 )
                [~,L] = lebesgue(x,d);   % lebesgue constant
                title(['Lebesgue constant = ' num2str(L)],FS,14)
            end
          end
          [gx,gy,button] = ginput(1);    % input new interpolation point
          if ( button == 3 ) break, end  % if right button, stop
          x = unique([x; gx]);           % #ok<AGROW>
        end
    end

%%
% P.S. If you want to see how a greedy person would choose their points,
% check out the Greedy Interpolation example [1].

end

%%
% References:
%
% [1] http://www.chebfun.org/examples/approx/html/GreedyInterp.shtml
