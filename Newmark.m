clc;clear;
% Newmark solver for a two degrees of freedom system under earthquake

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

M_mat=[m0 0; 0 m0];
K_mat=[2*k0 -k0; -k0 k0];

alpha1=2*w1*w2*(epsilon_1*w2-epsilon_2*w1)/(w2^2-w1^2);
alpha2=2*(epsilon_2*w2-epsilon_1*w1)/(w2^2-w1^2);

C_mat=alpha1*M_mat + alpha2*K_mat;

%%
% Newmark parameters
delta=0.50;
alpha=0.25*(0.5+delta)^2;

a0=1/(alpha*dt^2);
a1=delta/(alpha*dt);
a2=1/(alpha*dt);
a3=1/(2*alpha)-1;
a4=delta/alpha-1;
a5=0.5*dt*(delta/alpha-2);
a6=dt*(1-delta);
a7=dt*delta;

acceleration=zeros(n,1);
velocity=zeros(n,1);
displacement=zeros(n,1);

u(n,nt)=0;
v(n,nt)=0;
a(n,nt)=0;

for i=2:nt
    
    detF=-M_mat*input_acc(i)*[1; 1];
    
    effectivestiffness=a0*M_mat+a1*C_mat+K_mat;

    effectiveforce=M_mat*(a0*displacement+a2*velocity+a3*acceleration)+...
                   C_mat*(a1*displacement+a4*velocity+a5*acceleration)+...
                   detF;
    
    temp01=effectivestiffness\effectiveforce;
% temp01 = bicgstabl(effectivestiffness, effectiveforce, 1e-18, 200);

    
    incrementalsolution=temp01-displacement;
    temp02=a0*incrementalsolution-a2*velocity-a3*acceleration;
    temp03=velocity+a6*acceleration+a7*temp02;
    displacement=temp01; 
    acceleration=temp02;
    velocity=temp03;
    
    u(:,i)=displacement;
    v(:,i)=velocity;
    a(:,i)=acceleration;
    
    
end


%%
subplot(3,1,1)
plot(u(1,:));

subplot(3,1,2)
plot(v(1,:));

subplot(3,1,3)
plot(a(1,:));


