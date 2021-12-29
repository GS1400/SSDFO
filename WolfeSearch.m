
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [alpha,falp,g,exitflag,info] = WolfeSearch(fun,xinit,xRows, ...
%                           xCols,p,finit,fpinit,alpinit,tune,info)
%
% computes a steplength alpha > 0 that reduces the value of the 
% function fun along the search direction p, starting at xinit, 
% so that it satisfies the Wolfe conditions
%
% f(alpha) <= f(0) + rho*f'(0) 
% abs(f(alpha)) <= -sigma*f'(0).
%
% Here f(alpha) := fun(xinit + alpha*p), rho < 1/2  and rho < sigma < 1. 
%
% exitflag =  1: steplength alpha for which f(alpha) < fminimum was found 
% exitflag =  0: acceptable steplength was found
% exitflag = -1: nfmax reached
% exitflag = -2: secmax reached
% exitflag = -3: target accuracy reached
% exitflag = -4: no acceptable point could be found
%
% References:
%
% 1. R. Fletcher, Practical Methods of Optimization, John Wiley & Sons,
% 1987, second edition, section 2.6.
%
% 2. M. Al-Baali and R. Fletcher, An Efficient Line Search for 
% Nonlinear Least Squares, Journal of Optimization Theory and 
% Applications, 1986, Volume 48, Number 3, pages 359-377.
%
% 3. WolfeSearch has been rewitten according to line search used
% by fminunc
%    https://www.mathworks.com/help/optim/ug/fminunc.html
% 

function [alpha,falp,g,exitflag,info] = ...
WolfeSearch(fun,xinit,xRows,xCols,p,finit,fpinit,alpinit,tune,info)
        
TolFun        = tune.TolFun;
DiffMinChange = tune.DiffMinChange;
DiffMaxChange = tune.DiffMaxChange;
TypicalX      = tune.TypicalX;
rho           = tune.rho;
sigma         = tune.sigma;
fminimum      = tune.fminimum;
   

% Find a bracket of acceptable points
[a,b,fa,fpa,fb,fpb,alpha,falp,g,exitflag,info] = goodBrack(fun,...
xinit,xRows,xCols,p,finit,fpinit, alpinit,rho,sigma,fminimum, ...
DiffMinChange,DiffMaxChange, TypicalX,info);


if exitflag == 2 
  % goodBrack found a bracket containing acceptable points;
  % now getAccept finds acceptable point within bracket
    [alpha,falp,g,exitflag,info] = getAccept(fun,xinit,xRows,...
     xCols,p,finit,fpinit,a,b,fa,fpa,fb,fpb,rho,sigma,TolFun,...
     DiffMinChange, DiffMaxChange,TypicalX,info); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
