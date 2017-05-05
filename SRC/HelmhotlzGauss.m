%Scientific Computing MECE 5397
%Maria Violeta Paez 
%Project A - Helmholtz Equation
%Helmholtz Equation using Gauss-Seidel Method 

clc
clear all 

%% 
%Ghost nodes, made as input in order to manipulate the code
n=input('Enter your value for n= ')

%Given values, constants

gamma=-1; ax=-pi; ay=-pi; by=pi; bx=pi;

%Creating vector with linespace function
x=linspace(ax,bx,n); y=linspace(ay,by,n);

%Boundary conditions
u=zeros(n);
%B.C. are vectorized to optimize the code
% u(:,1)=ax; %(You can find it as a code in line 76)
u(:,n)=((bx-ax).^2.*cos((pi.*bx)./ax))+((y(:)-ay)./(by-ay)).*(bx.*(bx-ax).^2-((bx-ax).^2.*cos((pi.*bx)./ax)));
u(1,:)=x(:).*(x(:)-ax).^2;
u(n,:)=(x(:)-ax).^2.*cos(pi.*x(:)./ax);

h=bx/n; %Step Size 
err=1;
frequency=10; 
iter=0
tic;  %Timer to evaluate the performance 

%% Checkpointing
% Sometimes files take a long time to run to completion. As a result, sometimes they crash due to a variety of reasons: power failure, walltime limit, scheduled shutdown, etc. 
% Checkpoint/Restarting has long been a common technique to tackle this issue. Checkpointing/Restarting essentially means saving data to disk periodically so that, if need be, 
% you can restart the job from the point at which your data was last saved. 


% Before the start of the iteration loop, "check-in" each variable
% that should be checkpointed in the event of restarting the job

matfile = 'PoissonEquationSolution.mat';     % mandatory; name of checkpoint mat-file
s = struct();                                % mandatory; create struct for checkpointing
s = chkin(s,{'iter'});                       % mandatory; iter is iteration loop index
s = chkin(s,{'frequency'});                  % mandatory; frequency is checkpointing period 
                                             % i.e., how often to perform a save

% continue until all variables are checked in. Note that you are only
% checking in the variables, they don't need to have been already defined

chkNames = fieldnames(s);    % the full list of variables to checkpoint
nNames = length(chkNames);   % number of variables in list

%% 
while max(max(err(:)))>=1e-6  %Tolerance 
    iter=iter+1;
    
    % If you want to test the restart script, use the function pause(1) to slow down the while loop.
    % This will slow down the while loop to 1 sec per iteration so that ctrl + C can used be to
    % "kill" the code to simulate a computer crash. From there, use the restart script to restart the loop.  
    %pause(.05) Only used to with small values of n, to play w/ restart code
    if mod(iter, frequency) == 0 % If statement, checkpoints periodically (determined by the frequency)
        chkpt                    % chkpt script performs checkpointing (save) every *frequency* iterations
        fprintf(1, ['Checkpointing frequency is every %2d iterations.' ...
          'Data updated at iteration %3d\n'], ...
          frequency, iter);      % Confirm after each checkpointing event 
    end
    
    
    uold=u;
   
for  j=2:n-1
    for i=2:n-1
        F(i,j)=sin(pi.*((x(i)-ax)/(bx-ax))).*cos((pi/2).*(2.*(((y(j)-ay)/(by-ay))+1)));
        %Discritization 
        u(i,j)= (1/((gamma*h^2)-4))*((h^2)*F(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)));
    end
     k=j; i=j; j=1;
    u(k,1)= (1/((gamma*h^2)-4))*((h^2)*F(i,j)-(u(i+1,j)+u(i+1,j)+u(i,j+1)+u(i,j+1))); %This is the boundary conditions.
end
unew=u;
err=abs((uold-unew)./unew);
 fprintf(1, 'Completed iteration %d\n', iter);
end
timedoc=toc

%% Plot
figure
contourf(u)
colorbar('location','eastoutside','fontSize',12);
xlabel('X Number of Nodes in X-direction','fontSize',12);
ylabel('Y Number of Nodes in Y-direction','fontSize',12);
title('Gauss Seidel for Helmhotlz')
figure
mesh(x,y,u)
xlabel('X Number of Nodes in X-direction','fontSize',12);
ylabel('Y Number of Nodes in Y-direction','fontSize',12);
zlabel('Position U','fontSize',12);
title('Gauss Seidel for Helmhotlz');
