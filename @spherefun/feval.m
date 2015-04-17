function y = feval( f, theta, lambda ) 
% Evaluate a spherefun. 

lambda = lambda/pi - .5;   % painful... 
theta = theta/pi;
y = feval( f.Cols, lambda ) * f.BlockDiag * feval(f.Rows, theta).'; 

end 