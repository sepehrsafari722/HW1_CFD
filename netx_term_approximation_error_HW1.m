clc 
clear
close all
%input Datas
k=input('k= ');
h=input('h= ');
C=input('C= ');
A=input('A= ');
P=input('P= ');
T_0=input('T_0= ');
T_inf=input('T_inf= ');
L=input('L= ');
ro=input('density= ');
...Number of pieces (J) created on the length L
g=zeros(1,100);
G=zeros(1,100);
for w=1:100
J=500+(2*w+1);;
dx=L/J;
...Calculation m
m=h*P/(k*A);
...Calculation of coefficients a, b, d
a=1/dx.^2;
b=(-2/dx^2)-m;
...Create a coefficient matrix in left hand
s=zeros(J-1);
for i=1:J-1
    for j=1:J-1
        if i==j
            if i==J-1
               s(i,j)=b+(a/(1+(h.*dx/k))); 
            else
            s(i,j)=b;
            end
        elseif  i==j+1
            s(i,j)=a;
        elseif i==j-1
            s(i,j)=a;
        else 
            s(i,j)=0;
        end
    end
end
%%Create matrices with known values in right hand
q=zeros(J-1,1);
for j=1:J-1
    if j==1
        q(j,1)=-m*T_inf-(a*T_0);
    elseif j==J-1
        q(j,1)=-a*(h*dx*T_inf/(k+(h*dx)))-(m*T_inf);
    else
       q(j,1)=-m*T_inf;
    end
end
...Calculation of temperature distribution
T=linsolve(s,+q);
g(1,w)=dx;  
G(1,w)=m*(T(0.5*(J+1),1)-T_inf)*(dx^2)/2;
end
figure(1)
loglog(g,abs(G),'.')
title('Next Term Approximation Error')
xlabel('dx')
ylabel('Error')
figure(2)
slope=gradient(log(G))./gradient(log(g));
plot(g,slope,'.')
title('Slope error diagram')
xlabel('dx')
ylabel('Slope')