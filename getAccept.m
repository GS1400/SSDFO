%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getAccept %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [alpha,falp,g,exitflag,info] = getAccept(fun,xinit,xRows,...
% xCols,p,finit,fpinit,a,b,fa,fpa,fb,fpb,rho,sigma,...
% TolFun,DiffMinChange,DiffMaxChange,TypicalX,info) 
%
% getAccept finds an acceptable point alpha within a given bracket [a,b] 
% containing acceptable points. Notice that nf counts the total number 
% of function evaluations including those of goodBrack

function [alpha,falp,g,exitflag,info] = getAccept(fun,xinit,xRows,...
xCols,p,finit,fpinit,a,b,fa,fpa,fb,fpb,rho,sigma,...
TolFun,DiffMinChange,DiffMaxChange,TypicalX,info) 

 
tau2 = min(0.1,sigma); tau3 = 0.5; tol = TolFun/1000;

alpha = []; falp = []; g = []; 
while 1

    
    % pick alpha in reduced bracket
    brcktEndpntA = a + tau2*(b - a); brcktEndpntB = b - tau3*(b - a);
    % find global minimizer in bracket [brcktEndpntA,brcktEndpntB] 
    
    alpha = goodStep(brcktEndpntA,brcktEndpntB,a,b,fa,fpa,fb,fpb);
    
                                                                           
    % evaluate f(alpha) and f'(alpha)
    
    falp = fun(reshape(xinit(:)+alpha*p(:),xRows,xCols));
    if isnan(falp), falp= inf; end
    
   
    info.nf = info.nf +1;      
   
    % check stopping tests
    sec=(cputime-info.initTime);
    
    if sec>info.secmax, exitflag = -2; info.done =1; end
    
    if info.nf>=info.nfmax, exitflag = -1; info.done = 1; end
    
    
    info.qf  = (falp-info.fbest)/(info.finit-info.fbest);
    if info.qf<=info.accf, exitflag = -3;  info.done = 1; end
  
    info.sec = sec;
    
    
   fpalp = FFDD(fun,xinit(:)+alpha*p(:),falp,p,...
            DiffMinChange,DiffMaxChange); 
   info.nf   = info.nf + 1;
   
   if info.done, 
       g  = FFD(fun,xinit(:)+alpha*p(:),...
       falp,DiffMinChange,DiffMaxChange,TypicalX);
       info.nf   = info.nf + xRows;
       g = g(:); 
       return; 
   end  
   
    % check if roundoff errors are stalling convergence.
    
    if abs((alpha-a)*fpa)<=tol, 
       exitflag = -4; 
       return; 
    end
    
    
    

    % update bracket
    aold = a; bold = b; faold = fa; fbold = fb; 
    fpaold = fpa; fpbold = fpb;
    ok = (falp > finit + alpha*rho*fpinit || falp >= fa);
    if ok
        a = aold; b = alpha; 
        fa = faold; fb = falp;
        fpa = fpaold; fpb = fpalp;
    else
        ok = (abs(fpalp) <= -sigma*fpinit);
        if ok, 
            exitflag = 0; 
             g  = FFD(fun,xinit(:)+alpha*p(:),...
                  falp,DiffMinChange,DiffMaxChange,TypicalX);
             info.nf   = info.nf + xRows;
             g = g(:); 
            return; 
        end
        a = alpha; fa = falp; fpa = fpalp;
        ok = ((b - a)*fpalp >= 0);
        if ok
            b = aold; fb = faold; fpb = fpaold;
        else
            b = bold; fb = fbold; fpb = fpbold;
        end
    end
    
    % check if roundoff errors are stalling convergence
    if abs(b-a) < eps, 
        exitflag = -4; 
        return; 
    end      
    
end % of while




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%