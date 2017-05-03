%Scientific Computing MECE 5397
%Maria Violeta Paez 
%Project A - Helmholtz Equation
%Helmholtz Equation using Over-Relaxation Method

clc
clear all 

%Ghost nodes, made as input in order to manipulate the code
n=input('Enter your value for n= ')

%Given values, constants
gamma=-1; ax=-pi; ay=-pi; by=pi; bx=pi;
b=1.5; %Betta, for over relaxation 

%Creating vector with linespace function
x=linspace(ax,bx,n); y=linspace(ay,by,n);

%Boundary conditions
u=zeros(n); %Intial guess for Guass-Seidel
            %Method is zero for all interior nodes
u(:,1)=ax;
u(:,n)=((bx-ax).^2.*cos((pi.*bx)./ax))+((y(:)-ay)./(by-ay)).*(bx.*(bx-ax).^2-((bx-ax).^2.*cos((pi.*bx)./ax)));
u(1,:)=x(:).*(x(:)-ax).^2;
u(n,:)=(x(:)-ax).^2.*cos(pi.*x(:)./ax);

h=bx/n; %Step Size 
it=1000;
err=1;

while max(max(err(:)))>=1e-6
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
end
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