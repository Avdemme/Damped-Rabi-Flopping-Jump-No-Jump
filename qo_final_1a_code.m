%first we create our variables
syms psi_f psi_i psi_n ;
syms t dt t_rec step fns;
syms dp p p_rec p_avg;
syms hbar gamma;
syms H U L;
syms eps nor;
syms psq psq_avg msq var sdev sdev_tot;

%then we initialize our wvaefunction
psi_i = [0;1];
psi_f = [0;1];
psi_n = [0;1];

%then we adjust step sizes
dt = 0.0001;
t = 4;
step = t/dt;
fns = 100;

%then we create matrices to stor our probabilities
% and also a vector to represent a time scale
p_rec=zeros(step,fns);
t_rec=zeros(step,1);

%we create our time axis
for n=1:step
    t_rec(n) = (n-1)*dt;
end

%we initialize our probability
dp = 0;
p = 0;

%we set values for our constants
hbar = 1;
gamma = 1;

%we create our hamiltonian matrix
H = [-1*1i*hbar*gamma*(1/2), -3*hbar*gamma; -3*hbar*gamma, 0];

%from there we get our time evolutino matrix approximation
U = eye(2)-1i*(1/hbar)*dt*H;

%then we create our jump operators
L = 2*sqrt(gamma)*[0,0;1,0];
%and we must also initialize our random number
eps = 0;
%and initialize our norm number
nor = 0;

%initialize our standard deviation variables
psq=0;
psq_avg=0;
msq=0;
var=0;
sdev=0;
sdev_tot=0;

%now we must create multiple wave functions
for m=1:fns
    
%within this we must create a single quantum trajectory
%and we must reset for each one
    psi_i = [0;1];
    psi_f = [0;1];
    psi_n = [0;1];
    dp = 0;
    p = 0;
    eps = 0;
    nor = 0;
    
for n=1:step
    %for this code for every step we need to generate a random number
    eps = rand;
    
    %now we must make it so that we start with the wavefunction that 
    %was created at the end of the last step
    psi_i=psi_n;
    %now we calculate our dp
    dp = gamma*dt*(abs([1,0]*psi_i))^(2);
    %now we compare our random number to dp and decide how to 
    %evolve psi
    if eps < dp
        psi_f = L*psi_i;
    elseif eps > dp
        psi_f = U*psi_i;
    end
       
    %now we must renormalize our final state
    nor = norm(psi_f);
    psi_n= psi_f/nor;
    %now we calculate the probability of finding the normalized state in
    %the excited state
    p = (abs([1,0]*psi_n))^(2);
    %now we store that probability in a matrix to use for later
    p_rec(n,m) = p;
    
end
    
end
%now we must find the average probability for each step over all
%wavefunctions
p_avg=(1/fns)*sum(p_rec,2);

%plot(t_rec,p_avg);

%title('reproduction of Molmer Fig.3');
%xlabel('time(unit:1/gamma)')
%ylabel('Probability of the excited state');


%now we must calvulate the standard deviation so we start by squaring all
%the probabilities
psq = p_rec.^2;
psq_avg = (1/fns)*sum(psq,2);
msq = p_avg.^2;
var = psq_avg - msq;
sdev = var.^(1/2);
sdev_tot = (1/sqrt(fns))*sdev;
%now we plot our results

errorbar(t_rec,p_avg,sdev_tot);
title('reproduction of Molmer Fig.3 with error bars');
xlabel('time(unit:1/gamma)')
ylabel('Probability of the excited state');

