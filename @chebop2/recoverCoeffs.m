function p = recoverCoeffs(L)
%RECOVERCOEFFS   Recover coefficient functions of a linear operator.
%   P = RECOVERCOEFFS(L) returns, for a linear operator L, a CHEBFUN
%   quasimatrix P such that
%         Lu = P(:,1)*u + P(:,2)*u' + P(:,3)*u" + ... P(:,M+1)*u^(M),
%   where M is the difforder of the operator. If L is not linear, an error is
%   thrown.
%
%   For a block operator L, i.e., one defining a system of equations
%         Lu = [L_{1,1} L_{1,2} ... L_{1,S}] [ u_1 ]
%              [L_{2,1} L_{1,2} ... L_{1,S}] [ u_2 ]
%              [  ...     ...   ...   ...  ] [ ... ]
%              [L_{R,1} L_{R,2} ... L_{R,S}] [ u_S ],
%   P will be the RxS cell array such that P{J,K} = RECOVERCOEFFS(L_{J,K}).
%
%   [P L] = RECOVERCOEFFS(L) returns also the linop L, which can be useful if
%   the input was a linear chebop.
%
% Example 1:
%  [L x] = chebop(@(x,u) 0.5*diff(u,2) - sin(x).*diff(u) + x.*u);
%  p = recoverCoeffs(L)
%  norm(p - [x -sin(x) 0.5])
%
% Example 2:
%  [L x] = chebop(@(x,u) diff(sin(x).*(diff(cos(x).*u))),[-pi pi]);
%  p = recoverCoeffs(L)
%  norm(p - [-sin(2*x) 1-3*sin(x).^2 sin(2*x)/2])
%
% Example 3:
%  L = chebop(@(x,u,v) [diff(u,2), 0.5*diff(v)+exp(x)]);
%  p = recoverCoeffs(L)
%  norm([p{:}] - [0 0 1 0 0 0 .5])

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Convert to linop if input is a chebop. (But don't overwrite input as it's
% more efficient to evaluate the chebop .op than the linearised .oparray!)
if isa(L, 'chebop')
    L2 = linop(L);
else
    L2 = L;
end

% Initialize:
s = size(L2);                    % Determine the size of the system,
m = L2.diffOrder;                % and the difforder.
x = chebfun('x', L2.domain,2);   % Construct linear function on the domain,
x0 = chebfun(0, L2.domain);      % and the zero function.
p = cell(s);                     % Initialise output.
p0 = L*repmat(x0, 1, s(2));      % Compute non-autonomous component.

% The main routine:
for hh = 1:s(2)                 % Loop over each of the dependent variables.
    x0l = repmat(x0,1,hh-1);    % Set dep vars to the left to zero.
    x0r = repmat(x0,1,s(2)-hh); % Set dep vars to the right to zero.
    p1 = L*[x0l 1+0*x x0r];     % Evaluate all equations for [0 ... 1 ... 0]
    p1 = p1 - p0;               % Subtract non-autonomous component.
    for ll = 1:s(1)             % Loop over equations and assign.
        p{ll,hh} = p1{ll};
    end
    xk = x;                            % Update indep var to x.
    for kk = 1:max(m(:,hh))            % Loop over each x^k.
        tmp = L*[x0l xk(:,kk) x0r]-p0; % Evaluate for u = [0 ... x^k ... 0].
        for ll = 1:s(1)                % Loop over each equation.
            if kk > m(ll,hh)           % No coeffs of this order here.
                continue
            end 
            p{ll,hh}(:,kk+1) = tmp{ll}; % Assign the ll-th equation.
            for jj = 1:kk               % Extract the lower-order terms.
                p{ll,hh}(:,kk+1) = p{ll,hh}(:,kk+1) - p{ll,hh}(:,kk+1-jj).*xk(:,jj);
                p{ll,hh}(:,kk+1) = simplify(p{ll,hh}(:,kk+1)); % Simplify.
            end
        end
        xk = [xk x.*xk(:,end)/(kk+1)]; % Update indep var to x^k/k!
    end
end

end
