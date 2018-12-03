clc
clear all
close all
load Twod.mat
N=8 ;
Nv=2768;
M=7;

Xa=[x ones(Nv,1)];
X=Xa';
T=t';
r=zeros(N+1,N+1);
c=zeros(N+1,M);
W=zeros(N+1,M);

for i = 1: Nv
    
    for a=1:N+1         
        for b=1:N+1
            r(a,b)= r(a,b)+(X(a,i).*X(b,i))/Nv;
        end
    end
    
    for z=1:N+1
        for  y=1:M
            c(z,y)= c(z,y)+(X(z,i).*T(y,i))/Nv;
        end
    end
    
end

Ett=0;

for k=1:M
    tt=t(:,k);
    for p=1:Nv
        Ett= Ett+ (tt(p)^2)/Nv;
    end
    Et(k)=Ett;
    Ett=0;
end

R=r
C=c
Et;
 

for k=1:M
   for l=1:1 
    Nit=100;
    B2=0;
 W=zeros(N+1,1);
 disp(' i       Eg          XI           B2')

    for i=1:Nit
        a=zeros(N+1,1);b=zeros(N+1,1);
        gp=zeros(9,1);
        Num=0; Den=0;
        for n=1:N+1
            for m=1:N+1
                a(n)=a(n)+W(m).*R(n,m);
            end
        end
        for n=1:N+1
            g(n)=(-2*C(n,k)+2*a(n));
        end
         EG=0;
            for n=1:N+1
                EG=EG+g(n).*g(n);
            end
            
        XI=0;
        for n=1:N+1
            XI=XI+g(n).*gp(n);
            gp(n)=g(n);
        end
        
       
      
        for n=1:N+1
            
            Num = Num -g(n).*(C(n,k)-a(n));
            for m=1:N+1
                b(n) =b(n) +g(m).*R(n,m);
            end
            
            Den = Den +g(n).*b(n) ;
        end
        
       B2=Num/Den;
        if Den==0
            break
        end
        
        for n=1:N+1
            
            W(n)=W(n)-B2*g(n);
            
        end
       
        fprintf(' %d   %f   %.4f      %f\n', i, EG, XI, B2);

    end
    WW(k,:)=W;
    G(k,:)=g;
      end
   
end
XI;

WW 
q=linsolve(r,c)';

disp('-----------------------------------------')
MSE=0;
for i=1:M
   
    EX(i)=Et(i);
    for n=1:N+1
        EX(i)=EX(i)- 2*WW(i,n).*C(n,i);
    end
    
    for n=1:N+1
        for m=1:N+1
            
            EX(i)= EX(i)+ WW(i,n)*WW(i,m)*R(n,m);
            
        end
    end
    
    MSE= MSE+EX(i);
    fprintf('Error at node %d = %f\n', i, EX(i));
end
MSE
