%Scientific Computing MECE 5397
%Maria Violeta Paez 
%Project A - Helmholtz Equation
%Helmholtz Equation using Gauss-Seidel Method 

clc
clear all 

%Ghost nodes, made as input in order to manipulate the code
n=input('Enter your value for n= ')

%Given values, constants

gamma=-1; ax=-pi; ay=-pi; by=pi; bx=pi;

%Creating vector with linespace function
x=linspace(ax,bx,n); y=linspace(ay,by,n);

%Boundary conditions
u=zeros(n);
u(:,1)=ax;
u(:,n)=((bx-ax).^2.*cos((pi.*bx)./ax))+((y(:)-ay)./(by-ay)).*(bx.*(bx-ax).^2-((bx-ax).^2.*cos((pi.*bx)./ax)));
u(1,:)=x(:).*(x(:)-ax).^2;
u(n,:)=(x(:)-ax).^2.*cos(pi.*x(:)./ax);

h=bx/n; %Step Size 
it=1000;

for k=1:it
for  j=2:n-1
    for i=2:n-1
        F(i,j)=sin(pi.*((x(i)-ax)/(bx-ax))).*cos((pi/2).*(2.*(((y(j)-ay)/(by-ay))+1)));
        %Discritization 
        u(i,j)= 1/(4).*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)+F(i,j).*h.^2);
    end 
    u(i,n)= 1/(4).*(u(i-1,j)+u(i-1,j)+u(i,j-1)+u(i,j+1)+F(i,j).*h.^2);
end
end
%Plot 

figure
contourf(u)
colorbar('location','eastoutside','fontSize',12);
xlabel('X Number of Nodes in X-direction','fontSize',12);
ylabel('Y Number of Nodes in Y-direction','fontSize',12);
title('Gauss Seidel for Helmhotlz')
figure
surf(x,y,u)
xlabel('X Number of Nodes in X-direction','fontSize',12);
ylabel('Y Number of Nodes in Y-direction','fontSize',12);
zlabel('Position U','fontSize',12);
title('Gauss Seidel for Helmhotlz');
