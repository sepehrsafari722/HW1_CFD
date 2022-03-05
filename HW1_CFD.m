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
J=600;
dx=L/J;
...Calculation m
m=h*P/(k*A);
...Calculation of coefficients a, b, d
a=1/dx.^2;
b=(-2/dx^2)-m;
...Create a coefficient matrix
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
...ploting
figure(1)
T_end=(T(J-1,1)+(h*dx*T_inf/k))/(1+(h*dx/k));
T=[T_0;T;T_end]';
plot((0:J),T,'.')
title('Temperature Distribution')
xlabel('[0,L] , #points on [0,L]=600')
ylabel('Temperature')
...Error analysis
x=linspace(0,L,J+1);
%%T_exact=dsolve('D2T_exact-m*(T_exact-T_inf)','T_exact(0)=T_0','-k*DT_exact(L)=h*(T_exact(L)-T_inf)','x')
T_exact=T_inf + (exp(-m^(1/2)*x)*(T_0*h*exp(L*m^(1/2)) - T_inf*h*exp(L*m^(1/2)) + T_0*k*m^(1/2)*exp(L*m^(1/2)) - T_inf*k*m^(1/2)*exp(L*m^(1/2))))/(h*exp(L*m^(1/2)) - h*exp(-L*m^(1/2)) + k*m^(1/2)*exp(L*m^(1/2)) + k*m^(1/2)*exp(-L*m^(1/2))) - (exp(m^(1/2)*x)*(T_0 - T_inf)*(h*exp(-L*m^(1/2)) - k*m^(1/2)*exp(-L*m^(1/2))))/(h*exp(L*m^(1/2)) - h*exp(-L*m^(1/2)) + k*m^(1/2)*exp(L*m^(1/2)) + k*m^(1/2)*exp(-L*m^(1/2)));
figure(2)
loglog(x,abs(T-T_exact),'.')
title('Error between analytical and numerical value')
xlabel('Length [0,L]')
ylabel('Error')