% Script to run modified Traub-Miles hh model from Brette et al 2007
% Written by Dr Robert Stewart for Stewart & Bair, 2009
clear all; close all; warning('off','all'); format long g;

%Select sims
disp('Starting Hodgkin-Huxley simulation. Select simulation type:')
disp('1 = 1 spike current injection');
disp('2 = 10 spike current injection');
sim = input('3 = Analysis\n');

%set initial conditions for resting state
EL = -65; Vt = -63; E_alpha_n = Vt+15; E_beta_n = Vt+10; E_alpha_m = Vt+13;
E_beta_m = Vt+40; E_alpha_h = Vt+17; E_beta_h = Vt+40;
f0 = zeros(1,4); v = EL; f0(1) = v;

alpha_n = (0.032*5) * ((E_alpha_n-v)/5)/(exp((E_alpha_n-v)/5) - 1);
beta_n = 0.5*exp((E_beta_n-v)/40);
f0(2) = alpha_n/(alpha_n+beta_n);

alpha_m = (0.32*4) * ((E_alpha_m-v)/4)/(exp((E_alpha_m-v)/4) - 1);
beta_m = (0.28*5) * ((v-E_beta_m)/5)/(exp((v-E_beta_m)/5) - 1);
f0(3) = alpha_m/(alpha_m+beta_m);

alpha_h = 0.128 .* exp((E_alpha_h-v)./18); 
beta_h = 4 ./ (exp((E_beta_h-v)./5)+1);
f0(4) = alpha_h/(alpha_h+beta_h);

f0 = [f0,zeros(1,4)]; %Run substitutions
f0(5) = 0.5*exp((E_beta_n-f0(1))/40); %beta_n = 0.5*exp((E_beta_n-v)/40);
f0(6) = 0.128.*exp((E_alpha_h-f0(1))./18); %alpha_h=0.128.*exp((E_alpha_h-v)./18); 
f0(7) = exp((E_beta_h-f0(1))./5); %exp((E_beta_h-v)./5);
f0(8) = f0(4)/(f0(7)+1); %d = h/(c+1);

%Wide tolerance range like PS injection current sims
tols = 1*10.^(-2:-1:-16);%Wide tolerance range (15 values)
n_tols = length(tols); dt_ps = 1/10;
dt_vals = [1/100,1/200,1/400,1/600,1/800,1/1000,1/2000,1/4000,1/6000,1/8000,1/10000,1/20000,1/40000,1/60000,1/80000];

fp = zeros(100,1); ip = int32(fp); %floating point and integer parameter arrays
ip(end) = int32(200); %Order limit for Parker-Sochacki method
n_neurons = 10; n_ex = n_neurons; %number of neurons, number excitatory
t_end = 1000; sim_type = 0; ps_only = 0; syn_seed = 1; in_seed = 1;
temp = int32([syn_seed,sim_type,t_end,in_seed,ps_only,n_neurons]);
ip(1:length(temp)) = temp; 
I_inj = 0; tol = tols(1); dt_rk = dt_vals(1);
temp = [f0,I_inj,tol,dt_rk,dt_ps];
fp(1:length(temp)) = temp;
if(sim == 1)
  tm_inj(1,tols,n_tols,dt_vals,ip,fp);
elseif(sim == 2) 
  tm_inj(10,tols,n_tols,dt_vals,ip,fp);
elseif(sim==3)
  tm_inj_analysis;
end