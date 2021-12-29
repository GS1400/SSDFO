

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% driverSSDFO.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:3
    for j=1:100
        fprintf('=')
    end
    fprintf('\n')
end

fprintf(['This file illustrates the use of SSDFO by showing how ',...
    'to minimize the function f(x):=||Ax-b||_p^e \n',...
    'for p=2, e=1 in dimension n=10 from a random starting point. ',...
    'The goal assumed in the example is to\nreduce the',...
    ' initial objective function value by a factor of 1e-4 with',...
    ' at most nfmax=500*n function\nevaluations to return when one',...
    ' of the two limiting conditions is first satisfied.\n\n']);
disp('===============================================================')

clear;


solverPath = input(['Insert the SSDFO path \n',...
    '>> solverPath='],'s');


if ~exist(solverPath, 'dir')
    disp('the directory does not exist')
    return
end

disp('===============================================================')


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create problem definition


% define problem parameters (to be adapted to your problem)
n=10; % dimension
p=2;  % Norm in objective function
e=1;  % Exponent in objective function function

% create random matrix and right hand side for objective function
% (specific to the model problem; replace by whatever data you
% need to provide to your objectiv function)
A=rand(n)-0.5; 
b=-sum(A,2);

% create objective function f(x)=||Ax-b||_p^e
fun=@(x) norm(A*x-b,p).^e; 
% To solve your own problem, simply replace in the previous line 
% the expression after @(x) by your expression or function call. 
% Parameters in this expression are passed by value, hence are 
% fixed during minimization.

% start and stop info
x      = 2*rand(n,1);  % starting point

disp('===============================================================')


fullinfo=1;  Tuning=1;

% problem definition complete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem with SSDFO

tic; % set clock


if fullinfo    % define stop and print criteria
               % (indefinite run if no stopping criterion is given)
    % pass stop and print criteria
    % (indefinite run if no stopping criterion is given)
    
     st = struct('secmax',180,'nfmax',500*n,'finit',fun(x),...
         'fbest',0.001*fun(x),'accf',0.0001,'prt',2)
     
else
    st = []; % budgets are chosen inside SSDFO
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  solve the problem with SSDFO
if Tuning % self-tuning and info
    % given are the defaults
    % only the deviating parameters need to be set!
    
    tune = struct('DiffMaxChange',0.1,...
        'DiffMinChange',1e-8,...
        'TypicalX',ones(n,1),...
        'TolFun',0,...
        'mem',20,...
        'eps1',1e-8,...
        'gammaf3',1e-12,...
        'gammaf2',2,...
        'Deltaangle',1e-13,...
        'deltaa',1e-8)
    
else
   tune = []; % full tuning inside SSDFO is used
end
   

 % call SSDFO 
[x,f,info] = SSDFO(fun,x,st,tune);

if ~isempty(info.error),
   error = info.error
else
    % display output
    if st.prt>=0, 
        disp('SSDFO completed silently'); 
        disp(' ');
        disp('display output');
        disp(' ');
        info               % progress report by SSDFO
        nfused=info.nf    % number of function evaluations used
        secused=cputime-info.initTime   % time used up to now
        x                               % best point found     
        f                               % function value at x   
        fbest=st.fbest                  % target (for comparison)
    end
end

for i=1:3
    for j=1:100
        fprintf('=')
    end
    fprintf('\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
