
clc;clear;
% Runge Kutta solver for a two degrees of freedom system under earthquake

n=2; % story number
m0=1; % story mass
k0=2000*m0; % story stiffness
epsilon_1=0.05; epsilon_2=0.07; % damping ratios
w1=4.0; w2=10.0;  % 1st, 2nd circular frequency

nt=2000; % total time steps
dt=0.02; % time step size

load('waveinput5021-0.1g.txt');

input_acc=waveinput5021_0_1g(:,2);

%%
% construct mass, stiffness, damping matrices

M_mat=[m0 0; 0 m0];
K_mat=[2*k0 -k0; -k0 k0];

alpha1=2*w1*w2*(epsilon_1*w2-epsilon_2*w1)/(w2^2-w1^2);
alpha2=2*(epsilon_2*w2-epsilon_1*w1)/(w2^2-w1^2);

C_mat=alpha1*M_mat + alpha2*K_mat;

M_inv=[1/m0 0; 0 1/m0];

%%
u(n,nt)=0;
v(n,nt)=0;
a(n,nt)=0;


for i=2:nt
    
    acc_vec_n = input_acc(i-1)*[1; 1];
    acc_vec_nplus1 = input_acc(i)*[1; 1];
    acc_vec_nplushalf = (acc_vec_n + acc_vec_nplus1)/2;
    
    K1=v(:,i-1);
    M1=M_inv*(-C_mat*v(:,i-1)-K_mat*u(:,i-1)-M_mat*acc_vec_n);
    
    K2=v(:,i-1) + dt/2*M1;
    M2=M_inv*(-C_mat*(v(:,i-1) + dt/2*M1)-K_mat*(u(:,i-1) + dt/2*K1)-M_mat*acc_vec_nplushalf);
    
    K3=v(:,i-1) + dt/2*M2;
    M3=M_inv*(-C_mat*(v(:,i-1) + dt/2*M2)-K_mat*(u(:,i-1) + dt/2*K2)-M_mat*acc_vec_nplushalf);
    
    K4=v(:,i-1) + dt*M3;
    M4=M_inv*(-C_mat*(v(:,i-1) + dt*M3)-K_mat*(u(:,i-1) + dt*K3)-M_mat*acc_vec_nplus1);
    
    u(:,i) = u(:,i-1) + dt/6*( K1 + 2*K2 + 2*K3 + K4 );
    v(:,i) = v(:,i-1) + dt/6*( M1 + 2*M2 + 2*M3 + M4 );
    a(:,i) = ( v(:,i)-v(:,i-1) )/dt;
    
   
    
end

%%
subplot(3,1,1)
plot(u(1,:));

subplot(3,1,2)
plot(v(1,:));

subplot(3,1,3)
plot(a(1,:));


