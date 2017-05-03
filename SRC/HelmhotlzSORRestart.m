%% Restarting
matfile = 'PoissonEquationSolution';   % should match that in test_checkpoint.m
load(matfile);        % retrieve data from matfile
iter1 = iter+1;       % iter is the last time test_checkpoint issued
                      % a save; we start computing on the next step
%% Iterative Looping- SOR
while max(max(err(:)))>=1e-6  %Tolerance 
    iter=iter+1;
    
    % If you want to test the restart script, use the function pause(1) to slow down the while loop.
    % This will slow down the while loop to 1 sec per iteration so that ctrl + C can used be to
    % "kill" the code to simulate a computer crash. From there, use the restart script to restart the loop.  
    pause(.05)
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
        %Discritization, utilzing betta for over-relaxation
        u(i,j)= 1.*b/(4).*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+F(i,j).*h.^2)+(1-b).*u(i,j);
    end 
end
unew=u;
err=abs((uold-unew)./unew);
 fprintf(1, 'Completed iteration %d\n', iter);
end

timedoc=toc


%Plot

figure
contourf(u)
colorbar('location','eastoutside','fontSize',12);
xlabel('X Number of Nodes in X-direction','fontSize',12);
ylabel('Y Number of Nodes in Y-direction','fontSize',12);
title('SOR for Helmhotlz')
figure
surf(x,y,u)
xlabel('X Number of Nodes in X-direction','fontSize',12);
ylabel('Y Number of Nodes in Y-direction','fontSize',12);
zlabel('Position U','fontSize',12);
title('SOR for Helmhotlz');