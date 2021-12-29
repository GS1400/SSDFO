
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SSDFO.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,f,info] = SSDFO(fun,x,st,tune);
%
% solve the unconstrained black box optimization problem 
%    min f(x) 
%  
% fun      % function handle for f(.)
% x        % starting point (must be specified)
% st       % structure with stop and print criteria
%          % (indefinite run if no stopping criterion is given)
%  .secmax       %   stop if sec>=secmax (default: inf)
%  .nfmax        %   stop if nf>=nfmax   (default: inf)
%  .ftarget      %   stop if f<=ftarget  (default: -inf)
%  .prt          %   printlevel (default: -1)
%                %   -1: nothing, 0: litte, >=1: more and more
% tune     % optional structure containing tuning parameters
%          %   for details see below
%
% x        % best point found 
% f        % function value at best point found 
% info     % structure with performance info
%          %   for details see below
% 


function [x,f,info] = SSDFO(fun,x,st,tune)

info.error='';
% check function handle
if isempty(fun)
    
  message = 'SSDFO needs the function handle fun to be defined';
  info.error= message; 
  return
elseif ~isa(fun,'function_handle')
  message = 'fun should be a function handle';
  info.error= message;
  return
end


% starting point
if isempty(x)
  message = 'starting point must be defined';
  info.error= message;
  return      
elseif ~isa(x,'numeric')
  message = 'x should be a numeric vector'; 
  info.error= message;
  return       
end

[xRows,xCols] = size(x); % Store original user-supplied shape

% compute the initail function value
f = fun(x); info.nf=1;

if isnan(f), f= inf; end

info.finit=f;

%%%%%%%%%%%%%%%%%%%%%
% tuning parameters:
%%%%%%%%%%%%%%%%%%%%%
if ~exist('tune'), tune=[]; end;


% parameters for finite difference
% minimal value for the finite difference step size
if ~isfield(tune,'DiffMaxChange'), tune.DiffMaxChange = 0.1; end;
% maximal value for the finite difference step size
if ~isfield(tune,'DiffMinChange'), tune.DiffMinChange = 1e-8; end;

if ~isfield(tune,'TypicalX'), tune.TypicalX = ones(xRows,1); end;

% line search parameters
if ~isfield(tune,'rho'), tune.rho = 1e-4; end;
if ~isfield(tune,'sigma'), tune.sigma = 0.9; end;
if ~isfield(tune,'TolFun'), tune.TolFun = 0; end;
if ~isfield(tune,'fminimum'), tune.fminimum = f-1e8*(1+abs(f)); end

% parameters for the subspace
% subspace dimension
if ~isfield(tune,'mem'), tune.mem = 20; end;
% tiny parameter for adjusting the initial df
if ~isfield(tune,'eps1'), tune.eps1 = 1e-8; end;

% parameters for reducing df
if ~isfield(tune,'gammaf1'), tune.gammaf1 = 0.5; end;
if ~isfield(tune,'gammaf3'), tune.gammaf3 = 1e-12; end;

% parameter for expanding df
if ~isfield(tune,'gammaf2'), tune.gammaf2 = 2; end;

% tiny parameters for angle condition
% reqularization angle
if ~isfield(tune,'Deltaangle'), tune.Deltaangle = 1e-13; end;

% tiny factor for starting the step
if ~isfield(tune,'deltaa'), tune.deltaa = 1e-8; end;

if ~exist('st'), st=[]; end;


if isfield(st,'prt'), info.prt = st.prt; 
else, info.prt = 2; 
end;
% stopping criteria
if isfield(st,'secmax'), info.secmax=st.secmax;
else, info.secmax=inf;
end;
if isfield(st,'nfmax'), info.nfmax=st.nfmax;
else, info.nfmax=500*length(x);
end;

if isfield(st,'fbest'), info.fbest=st.fbest;
else, info.fbest=1e-4*fun(x);
end;

if isfield(st,'accf'), info.accf=st.accf;
else, info.accf = 1e-4;
end;

prt = info.prt;

info.initTime=cputime; info.done =0;

DiffMinChange = tune.DiffMinChange;
DiffMaxChange = tune.DiffMaxChange;
TypicalX      = tune.TypicalX;
deltaa        = tune.deltaa;


iter = 0; alpha = []; 


% Compute finite difference gradient at initial point, if needed
gradFd  = FFD(fun,x,f,DiffMinChange,DiffMaxChange,TypicalX);
info.nf = info.nf + xRows;

% the initial gradient vector
g = gradFd;
   
% Norm of initial gradient, used in stopping tests
normg = norm(g,Inf); 

% initialization for the subspace information
mem = min(tune.mem,xRows);
S   = zeros(xRows,mem);
Y   = zeros(xRows,mem);
H   = zeros(mem);
im  = 0; 
nh  = 0;

% initialize the initial value for df and dg
df  = tune.eps1*abs(f);

%%%%%%%%%%%%                    
% main loop
%%%%%%%%%%%%

done=0; info.deltaf=1; fbest=f;

while ~done
    
    iter = iter + 1;
    
     if prt>=0
       disp('=========================================================')
       disp(['Iteration ',num2str(iter)])
       disp('')
     end
    
    % compute the subspace direction
    if im == 0
       p = -g ; 
       if norm(p)==0
           p=rand(xRows,1)-0.5; gTp=-p'*p;
       else
          gTp=g'*p; 
       end
    else % subspace is not possible
     [p,gTp] = subspaceDir(S,Y,H,g,mem,nh,df,tune,info);
    end
     
    % perform line search along the subspace direction p
    alpha1 = 1;
    if iter == 1, alpha1 = min(1/normg,1); end  
    
    fold = f; gold = g;xold=x;

    % during line search, don't exceed the overall total maxFunEvals.
    info.maxFunEvals = info.nfmax - info.nf;
    
    WSok = (gTp < 0);
    if WSok
       [alpha,f,g,exitflag,info] = ... 
       WolfeSearch(fun,x,xRows,xCols,p,f,gTp,alpha1,tune,info);
   
   
       if info.done, break; end
       
      % if line search didn't finish successfully, find a heuristic
      % step size
      
        if exitflag == -4, 
            if prt>=0
                disp('line search did not finish successfully')
                disp('hence produce a heuristic step size')
            end
            WSok=0;
        end
    else
       if prt>=0
           disp('roundoff errors may be affecting convergence.')
           disp('hence produce a heuristic step size')
       end
    end
   
    
   if ~WSok
    
     x = xold; f = fold; 


        % find minimal step size
        if any(x==0 & p~=0)
          amin = deltaa*abs(f/gTp);
        else 
          ind=(p~=0);
          if any(ind)
            amin = deltaa*max(abs(f/gTp),min(abs(x(ind)./p(ind))));
          else  
            % zero direction found
            amin = 1;
          end
        end

       % compute the new point
       x  = x+amin*p;
       f  = fun(x); 
       if isnan(f), f= inf; end
       g = FFD(fun,x,f,DiffMinChange,DiffMaxChange,TypicalX);
       info.nf=info.nf+xRows+1;
       
      % check stopping test
      
      sec=(cputime-info.initTime);
      info.done=(sec>info.secmax)|(info.nf>=info.nfmax);
      info.qf  = (f-info.fbest)/(info.finit-info.fbest);
      info.done=info.done|(info.qf<=info.accf);
      info.sec = sec;
      
      
      
     if sec>info.secmax, exitflag = -2; end
    
     if info.nf>=info.nfmax, exitflag = -1; end
    
     if info.qf<=info.accf, exitflag = -3;  end
      
     if info.done, break; end
          
      s = amin*p;
   else
       % Update iterate
        s       = alpha*p;
        x       = x + s; 
        if prt>=0
            disp('line search finished successfully')
        end
   end
    
    % update df and dg
    
    if f<fold-df
       df = tune.gammaf1*(f-fold); 
    else
     df=max(tune.gammaf2*df,tune.gammaf3*(abs(f)+abs(fold)));
    end
    y        = g-gold;
    % update the best point
    if f<fbest, 
        
        fbest=f; 
        
        if prt>=0
          disp(['the best point was update; f = ',num2str(fbest)])
       end
    
    end;
    % Update the subspace information 
    [S,Y,H,nh,im] = updateSY(S,Y,H,nh,im,mem,s,y);
end % of while

x = reshape(x,xRows,xCols); % restore user shape
info.iterations = iter;
info.stepsize = alpha;
info.firstorderopt = norm(g,inf);

% determine status of SSDFO
switch exitflag
    case -3
       info.status = 'accuracy reached';
    case -1
       info.status = 'nfmax reached'; 
    case -2
       info.status = 'secmax reached';  
    otherwise 
       
end
info.status

