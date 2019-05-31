% Script to run benchmark Izhikevich model using Parker-Sochacki (PS),
% 4th order Runge-Kutta (RK) and Bulirsch-Stoer (BS) methods in iz_ps.c
% Written by Dr Robert Stewart for Stewart & Bair, 2009
clear all; close all;
warning('off','all'); format long;

%Select sims
disp('Starting Izhikevich simulation. Select simulation type:')
disp('1 = 1 spike current injection');
disp('2 = 10 spike current injection');
disp('3 = Current injection analysis');
disp('4 = Synaptic interations');
sim = input('5 = Synaptic interation analysis\n');

%Constant/default parameters
fp = zeros(100,1); ip = int32(fp); %floating point and integer parameter arrays
ip(end) = int32(200); %Order limit for Parker-Sochacki method

%Izhikevich version of the Benchmark neuron from Brette et al., 2007
C=200; vr=-65; vt=-50; k=1.3; a=0.03; b=-9.5; v_reset=-85; u_step=0; v_peak=48;
tau_ampa = 5; tau_gaba = 10; %time constants in ms
E_ampa = 0; E_gaba = -80; %Reversal potentials in mV
w_ampa = 6; w_gaba = -67; %Synaptic weights in nS
p_connect = 2.0/100.0; %Brette et al 2007 p391

tols = 1*10.^(-2:-1:-16);%Wide tolerance range (15 values)
n_tols = length(tols);
dt_vals = [1/4,1/6,1/8,1/10,1/20,1/40,1/60,1/80,1/100,1/200,1/400,1/600,1/800,1/1000,1/2000];

I_inj = 0; tol = 0; dt_rk = dt_vals(1); dt_ps = dt_vals(1); rand_inj_max = 200;
temp = [I_inj,tol,dt_rk,dt_ps,C,vr,vt,k,a,b,v_reset,u_step,v_peak,tau_ampa,tau_gaba,E_ampa,E_gaba,w_ampa,w_gaba,rand_inj_max,p_connect];
fp(1:length(temp)) = temp;

n_neurons = 1000; n_ex = n_neurons; t_end = 1000; 
sim_type = 0; syn_seed = 1; in_seed = 1; cnd = 1; ps_only = 0; trace_rec = 0;
temp = [n_neurons,sim_type,t_end,syn_seed,in_seed,cnd-1,ps_only,n_ex,trace_rec];
ip(1:length(temp)) = int32(temp);

if(sim == 1)
  iz_inj(1,tols,n_tols,dt_vals,ip,fp); %Run iz_inj for 1 spike
elseif(sim == 2) 
  iz_inj(10,tols,n_tols,dt_vals,ip,fp); %Run iz_inj for 10 spikes
elseif(sim == 3)
  iz_inj_analysis;
elseif(sim==4)
  neuron_vals = [1000,2000,4000]; %Set of network sizes
  iz_fix(neuron_vals,tols,dt_vals,ip,fp); %Full benchmarking simulation sweep, with multiple configs, inputs and conditions
elseif(sim==5)
  iz_fix_analysis(fp);  
end