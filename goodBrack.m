%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% goodBrack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [a,b,fa,fpa,fb,fpb,alpha,falp,g,exitflag,info] = goodBrack ...
% (fun,xinit,xRows,xCols,p,finit,fpinit,alpinit,rho,sigma,fminimum, ...
% DiffMinChange,DiffMaxChange,TypicalX,info) 
% 
% goodBrack finds a bracket [a,b] that contains acceptable points; 
% a bracket is the same as a closed interval, except that a > b is
% allowed.

function [a,b,fa,fpa,fb,fpb,alpha,falp,g,exitflag,info] = goodBrack ...
(fun,xinit,xRows,xCols,p,finit,fpinit,alpinit,rho,sigma,fminimum, ...
DiffMinChange,DiffMaxChange,TypicalX,info) 

tau1 = 9; % factor to expand the current bracket
a = []; b = []; fa = []; fpa = []; fb = []; fpb = []; g = []; 

% falp will contain f(alpha) for all trial points alpha
falp = finit;              

% fpalp will contain f'(alpha) for all trial points alpha
fpalp = fpinit;    

% set maximum value of alpha (determined by fminimum)
alpmax = (fminimum - finit)/(rho*fpinit); 

alpold = 0;

% first trial alpha is user-supplied
alpha = alpinit;
while 1

  % evaluate f(alpha) and f'(alpha)
  
  flod = falp; fpold = fpalp;
  
  
  falp = fun(reshape(xinit(:)+alpha*p(:),xRows,xCols));
  if isnan(falp), falp= inf; end
 
   info.nf   = info.nf+1;
  
    % check stopping tests
    sec=(cputime-info.initTime);
    
    if sec>info.secmax, exitflag = -2; info.done =1; end
    
    if info.nf>=info.nfmax, exitflag = -1; info.done = 1; end
    
    info.qf  = (falp-info.fbest)/(info.finit-info.fbest);
    
    if info.qf<=info.accf, exitflag = -3;  info.done = 1; end
  
    info.sec = sec;

   
   fpalp = FFDD(fun,xinit(:)+alpha*p(:),falp,p,...
                    DiffMinChange,DiffMaxChange); 
   info.nf  = info.nf + 1;
   if info.done, 
       g  = FFD(fun,xinit(:)+alpha*p(:),...
            falp,DiffMinChange,DiffMaxChange,TypicalX);
       info.nf   = info.nf + xRows;
       g = g(:); 
       return; 
   end  

   % terminate if f < fminimum
   if falp <= fminimum
      exitflag = 1;
       g  = FFD(fun,xinit(:)+alpha*p(:),...
            falp,DiffMinChange,DiffMaxChange,TypicalX);
       info.nf   = info.nf + xRows;
       g = g(:); 
     return 
   end
 
  % bracket located - case 1
  if falp > finit + alpha*rho*fpinit || falp >= flod    
    a = alpold; b = alpha;
    fa = flod; fpa = fpold;
    fb = falp; fpb = fpalp;
    exitflag = 2;
    g = []; 
    return 
  end
  
  % acceptable steplength found; no need to call sectioning phase
  if abs(fpalp) <= -sigma*fpinit
    exitflag = 0;
    g  = FFD(fun,xinit(:)+alpha*p(:),...
         falp,DiffMinChange,DiffMaxChange,TypicalX);
   info.nf   = info.nf + xRows;
   g = g(:); 
    return
  end
  
  % bracket located - case 2  
  if fpalp >= 0
    a = alpha; b = alpold;
    fa = falp; fpa = fpalp;
    fb = flod; fpb = fpold;
    exitflag = 2;
    g=[];
    return
  end
  
  % update alpha
  if 2*alpha - alpold < alpmax % if alpha + (alpha - alpold) < alpmax
      % brcktEndpntA = alpha + (alpha - alpold) >= alpmax
      brcktEndpntA = 2*alpha-alpold; 
      brcktEndpntB = min(alpmax,alpha+tau1*(alpha-alpold));
      % find global minimizer in bracket [brcktEndpntA,brcktEndpntB]
      % of 3rd-degree polynomial that interpolates f() and f'() 
      % at alpold and at alpha
      alpnew = goodStep(brcktEndpntA,brcktEndpntB,alpold,alpha,flod, ...
                                         fpold,falp,fpalp);
      alpold = alpha;
      alpha = alpnew;
  else
    alpha = alpmax;
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
