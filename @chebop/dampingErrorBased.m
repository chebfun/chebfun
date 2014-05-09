function [u, dampingInfo] = dampingErrorBased(N, u, rhs, delta, L, disc, dampingInfo)
%DAMPINGERRORBASED   Find the next step for damped Newton method.
%   DAMPINGERRORBASED finds the step-size lambda used in the damped Newton
%   iteration. It is affine invariant, and error controlled, that is, it seeks
%   to minimize the error in the SOLUTION SPACE, not the RESIDUAL SPACE.
%
%   [V, DAMPINGINFO] = DAMPINGERRORBASED(N, U, RHS, DELTA, L, DISC, DAMPINGINFO)
%   where
%    N:      Nonlinear CHEBOP
%    U:      Current guess of the solution of the BVP specified by N
%    RHS:    Current right-hand side of the differential equation
%    DELTA:  Current Newton corrections
%    L:      A LINOP, that is the linearization of N around U
%    DISC:   The CHEBDISCRETIZATION object arising from L
%    RHS:    Right hand side of ODE
%    V:      The new solution
%
%   Furthermore, the method takes in and returns the MATLAB struct DAMPINGINFO.
%   The fields of the struct are as follows:
%
%    errTol =       Description.
%    lambda =       Description.
%    lambdaMin =    Description. 
%    newtonCounter: Description.
%    normDelta:     Description.
%    normDeltaBar:  Description.
%    normDeltaOld:  Description.
%    deltaBar:      Description.
%   
%   For further details, see
%    [1] P. Deuflhard. Newton Methods for Nonlinear Problems. Springer, 2004.
%
%    [2] A. Birkisson. Numerical Solution of Nonlinear Boundary Value
%        Problems for Ordinary Differential Equations in the Continuous
%        Framework. DPhil Thesis, Oxford, 2013.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Describe the fields

% Extract info from the dampingInfo struct
errTol =        dampingInfo.errTol;
lambda =        dampingInfo.lambda;
lambdaMin =     dampingInfo.lambdaMin;
newtonCounter = dampingInfo.newtonCounter;
normDelta =     dampingInfo.normDelta;
normDeltaBar =  dampingInfo.normDeltaBar;
normDeltaOld =  dampingInfo.normDeltaOld;
deltaBar =      dampingInfo.deltaBar;

% Hard code x into N.op so that we don't have to keep recomputing it:
x = dampingInfo.x;
N.op = @(varargin) N.op(x, varargin{:});

% Monitors whether we want to accept the current steplength:
accept = 0;

% Indicates whether we are in prediction or correction mode (i.e,. whether we
% are finding the first value of lambda at a given Newton step, or whether we
% are correcting the value initially predicted for that step).
initPrediction = 1;

% Iterate until we find a step-size lambda that we accept:
while ( ~accept )
    
    % Check whether we want to predict a value for LAMBDA. Can only do so once
    % we have taken one Newton step, as it is based on information obtained from
    % the previous step
    if ( newtonCounter > 0 && initPrediction )
        % Compute a prediction value
        mu = (normDeltaOld*normDeltaBar) / ...
            (norm(deltaBar - delta, 'fro')*normDelta)*lambda;
        lambda = min(1, mu);
        % Indicate we will now be in correction mode until next Newton step.
        initPrediction = 0;
    end
    
    if ( lambda < lambdaMin )
        disp('Convergence failure')
        % Take full Newton step:
        uTrial = u + delta;
        accept = 1;
        normDeltaOld = normDelta;
        initPrediction = 1;
        lambda = 1;
        giveUp = 1;
        cFactor = NaN;
        continue
    end
    
    % Take a trial step
    uTrial = u + lambda*delta;
    % Evaluate the operator:
    NopTrial = feval(N, uTrial);
    
    % If N.op was stacked horizontally (i.e., V4 syntax), then we need to
    % convert NopTrial to a CHEBMATRIX:
    if ( isa(NopTrial, 'chebfun') && size(NopTrial, 2) > 1 )
        NopTrial = chebmatrix(cheb2cell(NopTrial).');
    end
    
    % Compute residual:
    deResFunTrial = NopTrial - rhs;
    
    % Compute a simplified Newton step using the current derivative of the
    % operator, but with a new right-hand side.
    [deltaBar, disc] = linsolve(L, deResFunTrial, disc);
    
    % We had two output arguments above, need to negate deltaBar:
    deltaBar = -deltaBar;    
    
    % TODO: Do we need to update the values of the RHS for the BCs here?
      
    % The norm of the simplified Newton step is used to compute a contraction
    % factor.
    normDeltaBar = norm(deltaBar);
    
    % Contraction factor:
    cFactor = normDeltaBar/normDelta;
    
    % Correction factor for the step-size:
    muPrime = (.5*normDelta*lambda^2) / norm(deltaBar-(1-lambda)*delta, 'fro');
    
    % If we don't observe contraction, decrease LAMBDA
    if ( cFactor >= 1 )
        lambda = min(muPrime, .5*lambda);
        % Go back to the start of the loop.
        continue
    end
    
    % New potential candidate for LAMBDA
    lambdaPrime = min(1, muPrime);
    
    if ( lambdaPrime == 1 && normDeltaBar < errTol )
        % TODO: We have converged within the damped phase! 
        % Do we need to treat this case separately within solvebvpNonlinear?
        u = uTrial + deltaBar; %#ok<NASGU>
        giveUp = 0; 
        break
    end
    
    % Switch to pure Newton if we are experiencing good convergence
    if ( lambdaPrime == 1 && cFactor < .5 )
        dampingInfo.damped = 0;
    end
    
    % TODO: Document
    if ( lambdaPrime >= 4*lambda )
        lambda = lambdaPrime;
        continue
    end
    
    % If we get all the way here, accept iterate, and tell the Newton iteration
    % to keep up the good work!
    accept = 1;
    giveUp = 0;
    
end

% Return UTRIAL:
u = uTrial;

% Update the dampingInfo structure:
dampingInfo.lambda =        lambda;
dampingInfo.cFactor =       cFactor;
dampingInfo.normDeltaBar =  normDeltaBar;
dampingInfo.deltaBar =      deltaBar;
dampingInfo.giveUp =        giveUp;

end
