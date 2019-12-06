clc
clear all
clear,
close


%%
%fp:passband edge 
%fa:stopband edge
%N: filter length 
%gam:gama variable
%b:Beta variable
%M:Number of frequency grids
%K:number of iteration for 
%Size of bounding box

%%

%Intializing variables
ri=ones(11,1)*0.005;
fp=0.4;
fa=0.5;
N=21;
M=100;
K=10000;
gam=.16;
b=0.08;
N1 = (N + 1)/2;

fp = fp*pi;
fa = fa*pi;
Mp = round(M*fp/(fp+pi-fa));
Ma = M - Mp;
f1= 0:fp/(Mp-1):fp;
f2 = fa:(pi-fa)/(Ma-1):pi;
f = [f1(:); f2(:)];

%Desired Omega values
omegad = [linspace(0,0.4*pi,44),linspace(0.5*pi,pi,56)];
omegad=f';


N=21;
K=(N-1)/2;
%Weight
W=ones(1,100);
%W(41:60)=0;
%intialize x0 with least-squares lowpass filter coeffcient with 0-0.4pi as  passband and 0.6pi to Pi stopband edges 
x0= fir1(20,0.45);


x = [x0(11) 2*x0((11+1):21)]';
x2=x;

%Desired Amplitude Ad
Ad=[ones(1,44),zeros(1,56)];

%Intialize C(w) and Ri*cos(w)
for i=1:length(omegad)
    
    for j=0:K
        c_omega(j+1,i)=cos(omegad(i)*j);
    end
    
    for j=0:K
        ricos(j+1,i)=ri(j+1)*abs( cos(j*omegad(i)) );
    end
    
end
   
 alpha(1)=0.1/(1+1);
 for k=1:1000000
     k
     
     %calculates E(x+d) for all values of Omega_Desired
     [omegastar(k),omegastar_ind(k)]=max( sum(ricos(:,:))'+abs(Ad' - c_omega(:,:)'*x(:,k))  );
     
     %Find the maximum delta point
     delta(:,k)=-sign(Ad(omegastar_ind(k)) - c_omega(:,omegastar_ind(k))'*x(:,k))*sign(c_omega(:,omegastar_ind(k))).*ri;
     
     y(k)= c_omega(:,omegastar_ind(k))'*x(:,k)  -  (Ad(omegastar_ind(k))  -  c_omega(:,omegastar_ind(k))'*delta(:,k));
     
     %Gradient function
     
     if y(k)>0
         grad_y(k)=1;
     elseif y(k)<0
         grad_y(k)=-1;
     elseif y(k)==0
         grad_y(k)= 0; 
     end
     
     
     g(:,k) = W(omegastar_ind(k))*c_omega(:,omegastar_ind(k))*(grad_y(k));
     
     %alpha beta Gamma calculation
     
     einfinty(k)= W(omegastar_ind(k)) * abs(c_omega(:,omegastar_ind(k))' * (x(:,k)+delta(:,k))  -  Ad(omegastar_ind(k)) );
     [minerror(k),minindex]=min(einfinty);
     
     beta(k)= 0.08/(k);
     gamma(k)=0.16/(k);
     alpha(k)=(einfinty(k)-minerror(k)+gamma(k))/((norm(g(:,k),2))^2);
     
     %Update the value of x
      if k==1
         x(:,k+1)= x(:,k)- alpha(k)*g(:,k);
      else
          x(:,k+1)= x(:,k)- alpha(k)*g(:,k) + beta(k)*(x(:,k) - x(:,k-1));
      end
     [minerror(k),minindex]=min(omegastar);
     
 end

 
xfinal=[ flipud(x(2:11,minindex))'/2 x(1,minindex)' (x(2:11,minindex)'/2) ];
 
omegad1 = [linspace(0,pi,1000)];
[h,w] = freqz(xfinal,1,1000);
plot(omegad1,(20*log10(abs(h))))
grid
title("Amplitude response of the robust minimax FIR filter")
ylabel("Amplitude")
xlabel("Normalized frequency")
 
figure
plot(20*log10(abs(c_omega(:,:)'*x(:,minindex))) )
grid
title("Amplitude response of the robust minimax FIR filter")
ylabel("Amplitude")
xlabel("Frequncy Grid")

iteration_num=[1:1000000]; 

figure
plot(iteration_num,einfinty)
axis([-100 10000 0 .4])
grid
title("Plot of error(norm infinity) with repesct to iteration")
ylabel("Error(norm infinity)")
xlabel(" iteration ")

iterationnum1=[1:200]
figure
plot(iterationnum1,einfinty(1:200))
grid
axis([-10 200 0 .4])
title("Plot of error(norm infinity) with repesct to first 200 iteration")
ylabel("error(norm infinity)")
xlabel(" iteration ")

figure
plot(xfinal)
grid
axis([0 22 -.2 0.6])
hold on
plot(x0)
legend('Intial coefficients of Filter', 'Final coefficients of Filter') 
title("Coefficients of Filter")
fprintf('The Coefficients of Filter are')
xfinal