function [x,fval,exitflag,output] = fMinSearch(funfcn,x,options,dis,varargin)
%FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
%   X = FMINSEARCH(FUN,X0) starts at X0 and attempts to find a local minimizer 
%   X of the function FUN.  FUN is a function handle.  FUN accepts input X and 
%   returns a scalar function value F evaluated at X. X0 can be a scalar, vector 
%   or matrix.
%
%   X = FMINSEARCH(FUN,X0,OPTIONS)  minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created
%   with the OPTIMSET function.  See OPTIMSET for details.  FMINSEARCH uses
%   these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, FunValCheck,
%   PlotFcns, and OutputFcn.
%
%   X = FMINSEARCH(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fminsearch' in PROBLEM.solver. 
%
%   [X,FVAL]= FMINSEARCH(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINSEARCH(...) returns an EXITFLAG that describes
%   the exit condition. Possible values of EXITFLAG and the corresponding
%   exit conditions are
%
%    1  Maximum coordinate difference between current best point and other
%       points in simplex is less than or equal to TolX, and corresponding 
%       difference in function values is less than or equal to TolFun.
%    0  Maximum number of function evaluations or iterations reached.
%   -1  Algorithm terminated by the output function.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINSEARCH(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name 
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fminsearch(@sin,3)
%     finds a minimum of the SIN function near 3.
%     In this case, SIN is a function that returns a scalar function value
%     SIN evaluated at X.
%
%     FUN can be an anonymous function:
%        X = fminsearch(@(x) norm(x),[1;2;3])
%     returns a point near the minimizer [0;0;0].
%
%     FUN can be a parameterized function. Use an anonymous function to
%     capture the problem-dependent parameters:
%        f = @(x,c) x(1).^2+c.*x(2).^2;  % The parameterized function.
%        c = 1.5;                        % The parameter.
%        X = fminsearch(@(x) f(x,c),[0.3;1])
%        
%   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.
%
%   See also OPTIMSET, FMINBND, FUNCTION_HANDLE.

%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
%   p.112-147, 1998.

%   Copyright 1984-2018 The MathWorks, Inc.
msg = 'ok';

defaultopt = struct('Display','notify','MaxIter','200*numberOfVariables',...
    'MaxFunEvals','200*numberOfVariables','TolX',1e-4,'TolFun',1e-4, ...
    'FunValCheck','off','OutputFcn',[],'PlotFcns',[]);

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && strcmpi(funfcn,'defaults')
    x = defaultopt;
    return
end

if nargin < 3, options = []; end

buildOutputStruct = nargout > 3;

% -------- Inputs are standardised for RAT, so no need to check ---------
% Detect problem structure input
% if nargin == 1
%     if isa(funfcn,'struct') 
%         [funfcn,x,options] = separateOptimStruct(funfcn);
%     else % Single input and non-structure
%         error('MATLAB:fminsearch:InputArg',...
%             sprintf('MATLAB:optimfun:fminsearch:InputArg'));
%     end
% end
% 
% if nargin == 0
%     error('MATLAB:fminsearch:NotEnoughInputs',...
%         sprintf('MATLAB:optimfun:fminsearch:NotEnoughInputs'));
% end
% 
% 
% % Check for non-double inputs
% if ~isa(x,'double')
%   error('MATLAB:fminsearch:NonDoubleInput',...
%     sprintf('MATLAB:optimfun:fminsearch:NonDoubleInput'));
% end
% -------------------------------------------------------------------

% n = numel(x);
% numberOfVariables = n;

% ------------- Check is done upstream ----------------
% Check that options is a struct
% if ~isempty(options) && ~isa(options,'struct')
%     error('MATLAB:fminsearch:ArgNotStruct',...
%         getString(message('MATLAB:optimfun:commonMessages:ArgNotStruct', 3)));
% end
% ------------------- AVH -----------------------------

% printtype = optimget(options,'Display',defaultopt,'fast');
tolx = optimget(options,'TolX',defaultopt,'fast');
tolf = optimget(options,'TolFun',defaultopt,'fast');
maxfun = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxiter = optimget(options,'MaxIter',defaultopt,'fast');
funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');

% In case the defaults were gathered from calling: optimset('fminsearch'):
% if ischar(maxfun) || isstring(maxfun)
%     if strcmpi(maxfun,'200*numberofvariables')
%         maxfun = 200*numberOfVariables;
%     else
%         error('MATLAB:fminsearch:OptMaxFunEvalsNotInteger',...
%             getString(message('MATLAB:optimfun:fminsearch:OptMaxFunEvalsNotInteger')));
%     end
% end
% if ischar(maxiter) || isstring(maxiter)
%     if strcmpi(maxiter,'200*numberofvariables')
%         maxiter = 200*numberOfVariables;
%     else
%         error('MATLAB:fminsearch:OptMaxIterNotInteger',...
%             getString(message('MATLAB:optimfun:fminsearch:OptMaxIterNotInteger')));
%     end
% end

switch dis      % Changed from TMW fminsearch
    case {'notify','notify-detailed'}
        prnt = 1;
    case {'none','off'}
        prnt = 0;
    case {'iter','iter-detailed'}
        prnt = 3;
    case {'final','final-detailed'}
        prnt = 2;
%     case 'simplex'
%         prnt = 4;
    otherwise
        prnt = 1;
end

% ----------------- Not using output functions for RAT ---------
% % Handle the output
% outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
% if isempty(outputfcn)
%     haveoutputfcn = false;
% else
%     haveoutputfcn = true;
%     xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
%     % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
%     outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
% end
% 
% % Handle the plot
% plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
% if isempty(plotfcns)
%     haveplotfcn = false;
% else
%     haveplotfcn = true;
%     xOutputfcn = x; % Last x passed to plotfcns; has the input x's shape
%     % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
%     plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
% end
% ---------------------- AVH ------------------------------

header = ' Iteration   Func-count     min f(x)         Procedure';

% Convert to function handle as needed.
% funfcn = fcnchk(funfcn,length(varargin));
% Add a wrapper function to check for Inf/NaN/complex values
controls = varargin{2};
problemStruct = varargin{1};
if funValCheck
    % Add a wrapper function, CHECKFUN, to check for NaN/complex values without
    % having to change the calls that look like this:
    % f = funfcn(x,varargin{:});
    % x is the first argument to CHECKFUN, then the user's function,
    % then the elements of varargin. To accomplish this we need to add the 
    % user's function to the beginning of varargin, and change funfcn to be
    % CHECKFUN.
    varargin = [{funfcn}, varargin];
    funfcn = @checkfun;
end
% 
n = numel(x);

% Initialize parameters
rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
onesn = ones(1,n);
% two2np1 = 2:n+1;
% one2n = 1:n;

% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = zeros(n,n+1); fv = zeros(1,n+1);
v(:,1) = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
x(:) = xin;    % Change x to the form expected by funfcn
[fv(:,1), result] = funfcn(x,varargin{:});
func_evals = 1;
itercount = 0;
coder.varsize('how',[1 Inf],[0 1]);
how = '';
% Initial simplex setup continues later

% Initialize the output and plot functions.
%
% ----------------------------------------
% RAT doesn't use output or plot functions...
%
% --------------------- AVH -----------

% if haveoutputfcn || haveplotfcn
%     [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'init',itercount, ...
%         func_evals, how, fv(:,1),varargin{:});
%     if stop
%         [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
%         if  prnt > 0
%             fprintf('%s \n', output.message)
%         end
%         return;
%     end
% end

% Print out initial f(x) as 0th iteration
if prnt == 3
    triggerEvent(coderEnums.eventTypes.Message, sprintf('\n%s\n', header));
    triggerEvent(coderEnums.eventTypes.Message, ...
                 sprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how));
% elseif prnt == 4
    % Option never used in RAT
    
    
%     formatsave.format = get(0,'format');
%     formatsave.formatspacing = get(0,'formatspacing');
%     % reset format when done
%     oc1 = onCleanup(@()set(0,'format',formatsave.format));
%     oc2 = onCleanup(@()set(0,'formatspacing',formatsave.formatspacing));
%     format compact
%     format short e
%     fprintf('%s \n', ' ')
%     fprintf('%s \n', how)
%     fprintf('%s \n', 'v = ')
%     fprintf('%g \n', v)
%     fprintf('%s \n', 'fv = ')
%     fprintf('%g \n', fv)
%     fprintf('%s \n', 'func_evals = ')
%     fprintf('%g \n', func_evals)
end

triggerEvent(coderEnums.eventTypes.Plot, result, problemStruct);

% OutputFcn and PlotFcns call
% if haveoutputfcn || haveplotfcn
%     [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
%         func_evals, how, fv(:,1),varargin{:});
%     if stop  % Stop per user request.
%         [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
%         if  prnt > 0
%             fprintf('%s \n', output.message)
%         end
%         return;
%     end
% end

% Continue setting up the initial simplex.
% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = 0.05;             % 5 percent deltas for non-zero terms
zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
for j = 1:n
    y = xin;
    if y(j) ~= 0
        y(j) = (1 + usual_delta)*y(j);
    else
        y(j) = zero_term_delta;
    end
    v(:,j+1) = y;
    x(:) = y; [f, result] = funfcn(x,varargin{:});
    fv(1,j+1) = f;
end

% sort so v(1,:) has the lowest function value
[fv,j] = sort(fv);
v = v(:,j);

how = 'initial simplex';
itercount = itercount + 1;
func_evals = n+1;
if prnt == 3 && rem(itercount, controls.updateFreq) == 0
    triggerEvent(coderEnums.eventTypes.Message, ...
                 sprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how));
% elseif prnt == 4
%     fprintf('%s \n', ' ')
%     fprintf('%s \n', how)
%     fprintf('%s \n', 'v = ')
%     fprintf('%g \n', v)
%     fprintf('%s \n', 'fv = ')
%     fprintf('%g \n', fv)
%     fprintf('%s \n', 'func_evals = ')
%     fprintf('%g \n', func_evals)
end
if rem(itercount, controls.updatePlotFreq) == 0
    triggerEvent(coderEnums.eventTypes.Plot, result, problemStruct);
end
if isRATStopped(controls.IPCFilePath)
    [x, fval, exitflag, output] = cleanUpInterrupt(v(:,1), fv(:,1), itercount, func_evals, prnt);
    return
end
% OutputFcn and PlotFcns call
% if haveoutputfcn || haveplotfcn
%     [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
%         func_evals, how, fv(:,1),varargin{:});
%     if stop  % Stop per user request.
%         [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
%         if  prnt > 0
%             fprintf('%s \n', output.message)
%         end
%         return;
%     end
% end
% exitflag = 1;

% Main algorithm: iterate until 
% (a) the maximum coordinate difference between the current best point and the 
% other points in the simplex is less than or equal to TolX. Specifically,
% until max(||v2-v1||,||v3-v1||,...,||v(n+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the 
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)
% The iteration stops if the maximum number of iterations or function evaluations 
% are exceeded
while func_evals < maxfun && itercount < maxiter
    if max(abs(fv(1)-fv(2:n+1))) <= max(tolf,10*eps(fv(1))) && ...
            max(max(abs(v(:,2:n+1)-v(:,onesn)))) <= max(tolx,10*eps(max(v(:,1))))
        break
    end
    
    % Compute the reflection point
    
    % xbar = average of the n (NOT n+1) best points
    xbar = sum(v(:,1:n), 2)/n;
    xr = (1 + rho)*xbar - rho*v(:,end);
    x(:) = xr; [fxr, result] = funfcn(x,varargin{:});
    func_evals = func_evals+1;
    
    if fxr < fv(:,1)
        % Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
        x(:) = xe; [fxe, result] = funfcn(x,varargin{:});
        func_evals = func_evals+1;
        if fxe < fxr
            v(:,end) = xe;
            fv(:,end) = fxe;
            how = 'expand';
        else
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        end
    else % fv(:,1) <= fxr
        if fxr < fv(:,n)
            v(:,end) = xr;
            fv(:,end) = fxr;
            how = 'reflect';
        else % fxr >= fv(:,n)
            % Perform contraction
            if fxr < fv(:,end)
                % Perform an outside contraction
                xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
                x(:) = xc; [fxc, result] = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                
                if fxc <= fxr
                    v(:,end) = xc;
                    fv(:,end) = fxc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xcc = (1-psi)*xbar + psi*v(:,end);
                x(:) = xcc; [fxcc, result] = funfcn(x,varargin{:});
                func_evals = func_evals+1;
                
                if fxcc < fv(:,end)
                    v(:,end) = xcc;
                    fv(:,end) = fxcc;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                for j=2:n+1
                    v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
                    x(:) = v(:,j); [fv(:,j), result] = funfcn(x,varargin{:});
                end
                func_evals = func_evals + n;
            end
        end
    end
    [fv,j] = sort(fv);
    v = v(:,j);
    itercount = itercount + 1;
    if prnt == 3 && rem(itercount, controls.updateFreq) == 0
        triggerEvent(coderEnums.eventTypes.Message, sprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how));
%     elseif prnt == 4
%         fprintf('%s \n', ' ')
%         fprintf('%s \n', num2str(how))
%         fprintf('%s \n', 'v = ')
%         fprintf('%s \n', v)
%         fprintf('%s \n', 'fv = ')
%         fprintf('%s \n', fv)
%         fprintf('%s \n', 'func_evals = ')
%         fprintf('%s \n', num2str(func_evals))
    end
    if rem(itercount, controls.updatePlotFreq) == 0   
        triggerEvent(coderEnums.eventTypes.Plot, result, problemStruct);
    end
    if isRATStopped(controls.IPCFilePath)
        [x, fval, exitflag, output] = cleanUpInterrupt(v(:,1), fv(:,1), itercount, func_evals, prnt);
        return
    end
    % OutputFcn and PlotFcns call
%     if haveoutputfcn || haveplotfcn
%         [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',itercount, ...
%             func_evals, how, fv(:,1),varargin{:});
%         if stop  % Stop per user request.
%             [x,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
%             if  prnt > 0
%                 fprintf('%s \n', output.message)
%             end
%             return;
%         end
%     end
end   % while

x(:) = v(:,1);
fval = fv(:,1);

if prnt == 3 && rem(itercount, controls.updateFreq) ~= 0
    % This should ensure the final result is printed at the end of a run irrespective of update frequency
    triggerEvent(coderEnums.eventTypes.Message, sprintf(' %5.0f        %5.0f     %12.6g         %s\n', itercount, func_evals, fv(1), how));
end
if rem(itercount, controls.updatePlotFreq) ~= 0
    % This should ensure the final result is always plotted irrespective of update frequency
   triggerEvent(coderEnums.eventTypes.Plot, result, problemStruct);
end

% OutputFcn and PlotFcns call
% if haveoutputfcn || haveplotfcn
%     callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'done',itercount, func_evals, how, fval, varargin{:});
% end

if func_evals >= maxfun
    printMsg = prnt > 0;
    if buildOutputStruct || printMsg
        %msg = getString(message('MATLAB:optimfun:fminsearch:ExitingMaxFunctionEvals', sprintf('%f',fval)));
        msg = sprintf('Exiting: Max function evals reached');
    end
    exitflag = 0;
elseif itercount >= maxiter
    printMsg = prnt > 0;
    if buildOutputStruct || printMsg
        %msg = getString(message('MATLAB:optimfun:fminsearch:ExitingMaxIterations', sprintf('%f',fval)));
        msg = sprintf('Exiting: Max iterations reached');
    end
    exitflag = 0;
else
    printMsg = prnt > 1;
    if buildOutputStruct || printMsg
        msg = sprintf('Exiting - X satisfies termination criteria: TolX %e, TolF %e',tolx,tolf);

    end
    exitflag = 1;
end

if buildOutputStruct
    output.iterations = itercount;
    output.funcCount = func_evals;
    output.algorithm = 'Nelder-Mead simplex direct search';
    output.message = msg;
end

if printMsg
    triggerEvent(coderEnums.eventTypes.Message, sprintf('\n%s\n', msg));
end
end
%--------------------------------------------------------------------------
% function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,state,iter,...
%     numf,how,f,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.
%
% state - can have the values 'init','iter', or 'done'.

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
% optimValues.iteration = iter;
% optimValues.funccount = numf;
% optimValues.fval = f;
% optimValues.procedure = how;

% xOutputfcn(:) = x;  % Set x to have user expected size
% stop = false;
% state = char(state);
% Call output functions

% ---- Remove these from function for compile - AVH
% if ~isempty(outputfcn)
%     switch state
%         case {'iter','init'}
%             stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
%         case 'done'
%             callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
%     end
% end
% % Call plot functions
% if ~isempty(plotfcns)
%     switch state
%         case {'iter','init'}
%             stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
%         case 'done'
%             callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
%     end
% end
% -----------------------------------

%--------------------------------------------------------------------------
function [x, fval, exitflag, output] = cleanUpInterrupt(optX, optVal, iteration, funccount, display)
    x = optX;
    fval = optVal;
    exitflag = -1;
    output.iterations = iteration;
    output.funcCount = funccount;
    output.algorithm = 'Nelder-Mead simplex direct search';
    output.message = 'Optimisation terminated by user';
    if display > 0
        triggerEvent(coderEnums.eventTypes.Message, sprintf('\n%s\n', output.message));
    end
end

%--------------------------------------------------------------------------
% function f = checkfun(x,userfcn,varargin)
% CHECKFUN checks for complex or NaN results from userfcn.

% f = userfcn(x,varargin{:});
% Note: we do not check for Inf as FMINSEARCH handles it naturally.
% if isnan(f)
%     error('MATLAB:fminsearch:checkfun:NaNFval','Target function is NaN');
% elseif ~isreal(f)
%     error('MATLAB:fminsearch:checkfun:ComplexFval',...
%         getString(message('MATLAB:optimfun:fminsearch:checkfun:ComplexFval', localChar( userfcn ))));  
%         error(sprintf('Target function is complex'));
% end

%--------------------------------------------------------------------------
% function strfcn = localChar(fcn)
% % Convert the fcn to a character array for printing
% 
% if ischar(fcn)
%     strfcn = fcn;
% elseif isstring(fcn) || isa(fcn,'inline')
%     strfcn = char(fcn);
% elseif isa(fcn,'function_handle')
%     strfcn = func2str(fcn);
% else
%     try
%         strfcn = char(fcn);
%     catch
%         strfcn = getString(message('MATLAB:optimfun:fminsearch:NameNotPrintable'));
%     end
% end


