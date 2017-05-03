%Maria Violeta Paez 
%Helmholtz Equation - MECE 5397
clc
clear all 
%Ghost nodes, made as input in order to manipulate the code
n=input('Enter your value for n')
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
h=length(x)/n

for  j=2:n-1
    for i=2:n-1
        F=sin(pi.*((x(i)-ax)/(bx-ax))).*cos((pi/2).*(2.*(((y(j)-ay)/(by-ay))+1)));
        u(i,j)= 1/(4-gamma)*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))*(1/h);
    end 
end
surf(x,y,u)