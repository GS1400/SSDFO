%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% subspaceDir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [p,gTp] = subspaceDir(S,Y,H,g,mem,nh,df,tune,info)
%
% compute new subspace direction
%

function [p,gTp] = subspaceDir(S,Y,H,g,mem,nh,df,tune,info)

if nh<mem, S=S(:,1:nh); Y=Y(:,1:nh); H=H(1:nh,1:nh); end

% compute Lam
Lam  = sqrt(sum(Y.^2'))./sqrt(sum(S.^2'));
 
% Adjust the components of Lam
Ok   = find(isnan(Lam) | Lam==0 | isinf(Lam));
Lam(Ok) = 1;  Lam = Lam(:); 

% find z by solving the linear system Hz=-c
c      = S'*g;
warning off
z   = -H\c;
warning on

gamma1 = c'*z;
gamma2 = 0.5*z'*H*z;
beta = -gamma1/gamma2;
              
% check whether beta is real or not
ok     = (isreal(beta)& ~isinf(beta)& ~isnan(beta));
ok=ok&gamma1<0&gamma2>0;
if ok 
   sub=1; beta = min(1,beta);
   if ~isnan(beta) & ~isinf(beta) & beta>0
       % compute the model function
        falp = beta*gamma1+0.5*beta^2*gamma2;
        % check whether both the function and gradient model
        % are decreased or not
        ok = (falp <= -df );
        ok = ok & isreal(falp)&isnan(falp); 
        if ok
            z   = beta*z; p =S*z;
            if any(isnan(p)|isinf(p))|norm(p)==0, sub=0; end
        else
            sub=0;
        end
   else
       sub=0;
   end
else
    sub=0;
end


if ~sub
    % compute the matrices U and M
    SS = S; on = ones(size(SS,2),1);

    U = Y- Lam(:,on).*SS; M = (Y'*(Y./Lam(:,on)))-H;

    % find z
    warning off
    z   = M\(U'*(g./Lam));
    warning on

    % chceck whether z is a zero vector or has been contaminated 
    % by NaN or inf 
    nok = (any(isnan(z)|isinf(z))|norm(z)==0);

    if nok, % z is not suitable;
            % new scaled steepest descent direction is constructed
       if info.prt>=1, 
           disp(['z was contaminated by nan or inf;',...
                 ' scaled steepest descent is tried']); 
       end
       p=-g./Lam;
    else 
          % limited memory direction is constructed
          if info.prt>=1, 
            disp('limited memory direction is constructed'); 
          end
          p = (U*z-g)./Lam; 
    end

end
    
gTp = g'*p;

% Moving away from maximazer or saddle point
if ~isnan(gTp)& gTp~=0
   if gTp>0, I = find(g.*p>0); p(I)=-p(I);  gTp = g'*p; end;
   % Enforce angle condition

    sigma1   = g'*g; sigma2  = p'*p;  sigma12  = sigma1*sigma2;
    sigmanew = gTp/sqrt(sigma12);

    angleOk = (sigmanew<=-tune.Deltaangle& isreal(sigmanew)&...
        ~isnan(sigmanew)&sigmanew<0);

    if angleOk, 
       % the angle condition holds; the direction does need
       % to be modeifed
       if info.prt>=1, 
           disp('the angle condition holds'); 
       end   
    else % the angle condition does not hold
       if info.prt>=1, 
           disp('enforce the angle condition'); 
       end   
       w    = (sigma12*max(eps,1-sigmanew^2))/(1-tune.Deltaangle^2);
       t    = (gTp+tune.Deltaangle*sqrt(w))/sigma1;
       wtOk = (w>0 & isfinite(t)&~isnan(t));
       if wtOk, p = p-t*g; gTp = g'*p;
           if gTp>0, I = find(g.*p>0); p(I)=-p(I); gTp = g'*p; end
       else 
           p = -g./Lam;  gTp = g'*p;
       end
    end 
else
    p=-g; gTp = g'*p;
end

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%