% Recurrent network simulation function, with fixed 1ms delays
% Written by Dr Robert Stewart for Stewart & Bair, 2009
function [] = iz_fix(neuron_vals,tols,dt_vals,ip,fp)
  file_stub = 'iz_fix';
  %Benchmarking simulation sweep, with multiple configs, inputs and conditions
  sim_type = 2; ip(2) = sim_type; %fixed 1ms synaptic event delays
  
  syn_seed = input('Select Network Config: 1-5 (All 5 needed before analysis)');
  if(syn_seed<1 || syn_seed>5)
    error('Network configuration must be between 1 and 5')
  end
  ip(4) = uint32(syn_seed);

  w_ampa_orig = fp(18); w_gaba_orig = fp(19);
  tols_orig = tols; dt_vals_orig = dt_vals; %store original values
  %Select representative tolerances from range used in current injection sims
	dt_vals = dt_vals([1,9,end,end]);
	tols = [tols([1,9,end]),0];
	n_cnd = 4;
	in_seeds = 1:10; %Input pattern indices
  conditions = [1,2,3,4]; %Four error tolerance conditions	
  t_end = ip(3); %sim duration in ms
	max_spikes = 200000;
  vr = fp(6);
  for n_neurons = neuron_vals
    if(n_neurons<4000)
      trace_rec=0;
    else
      trace_rec=1;
    end
  	fac = 4000/n_neurons;
  	weight_gain = sqrt(fac);
  	weight_gain=1.4; %Gain for roughly 10 spikes/cell 
  	w_gaba = w_gaba_orig/weight_gain;
  	w_ampa = w_ampa_orig/weight_gain;	
	  n_ex = n_neurons*.8; %number of neurons, number excitatory
    
    %Put parameters into fp
    ip(1) = uint32(n_neurons); ip(8) = uint32(n_ex); ip(9) = trace_rec;
    fp(18) = w_ampa; fp(19) = w_gaba;
	  
	  all_RK_tf = zeros(max_spikes,n_cnd);
	  all_RK_nrn = zeros(max_spikes,n_cnd);
	  all_PS_tf = zeros(max_spikes,n_cnd);
	  all_PS_nrn = zeros(max_spikes,n_cnd);
	  all_BS_tf = zeros(max_spikes,n_cnd);
	  all_BS_nrn = zeros(max_spikes,n_cnd);
	  all_t_cpu = zeros(3,n_cnd);
	  all_n_spikes = zeros(3,n_cnd);
	  all_RK_v = zeros(t_end,n_cnd);
	  all_PS_v = zeros(t_end,n_cnd);
	  all_BS_v = zeros(t_end,n_cnd);
	  disp(['Config seed = ',num2str(syn_seed)])
	  for in_seed = in_seeds
	    disp(['Input seed = ',num2str(in_seed)])
      ip(5) = int32(in_seed); 
	    for cnd = conditions
	    	disp(['Condition = ',num2str(cnd)])
        ip(6) = int32(cnd-1);
	    	if(cnd==4) ps_only = 1; %Run PS alone at condition 4 (reference solution)
	    	else ps_only = 0;
        end
        ip(7) = int32(ps_only); 
	    	
	  		%Set error tolerance and RK dt given condition
	  		tol = tols(cnd); dt_rk = dt_vals(cnd);
        fp(2:3) = [tol,dt_rk]; 
	      [RK_tf,RK_nrn,PS_tf,PS_nrn,BS_tf,BS_nrn,t_cpu,RK_v,PS_v,BS_v,i_stats,f_stats] = iz_ps(fp,ip);
	      ps_order_mean(cnd) = f_stats(1)
	      sub_mean(cnd) = f_stats(2)
	      ps_order_max(cnd) = i_stats(1)
	      sub_max(cnd) = i_stats(2)
	      bs_fails(cnd) = i_stats(3)
	            
	      RK_ind = find(RK_tf); RK_tf = RK_tf(RK_ind); RK_nrn = RK_nrn(RK_ind);
	      BS_ind = find(BS_tf); BS_tf = BS_tf(BS_ind); BS_nrn = BS_nrn(BS_ind);
	      PS_ind = find(PS_tf); PS_tf = PS_tf(PS_ind); PS_nrn = PS_nrn(PS_ind);
	      RK_n_spikes = length(RK_ind);
	      BS_n_spikes = length(BS_ind);
	      PS_n_spikes = length(PS_ind);
	      n_spikes = [RK_n_spikes,BS_n_spikes,PS_n_spikes];
	      disp(['Total number of spikes [RK BS PS]: ',num2str(n_spikes)])  
	      %mean(n_spikes) 
	      
	      %sort spikes into strict temporal order
	      [RK_tf,ind] = sort(RK_tf);
	      RK_nrn = RK_nrn(ind);
	      [BS_tf,ind] = sort(BS_tf);
	      BS_nrn = BS_nrn(ind);
	      [PS_tf,ind] = sort(PS_tf); 
	      PS_nrn = PS_nrn(ind);
	      
	      %gather results across conditions
	      all_RK_tf(1:n_spikes(1),cnd) = RK_tf;
	      all_RK_nrn(1:n_spikes(1),cnd) = RK_nrn;
	      all_BS_tf(1:n_spikes(2),cnd) = BS_tf;
	      all_BS_nrn(1:n_spikes(2),cnd) = BS_nrn;  
	      all_PS_tf(1:n_spikes(3),cnd) = PS_tf;
	      all_PS_nrn(1:n_spikes(3),cnd) = PS_nrn;
	      
	      all_t_cpu(:,cnd) = t_cpu;
	      all_n_spikes(:,cnd) = n_spikes';
	      
	      all_RK_v(:,cnd) = RK_v+vr;
	      all_BS_v(:,cnd) = BS_v+vr;
	      all_PS_v(:,cnd) = PS_v+vr;
	    end
	    all_t_cpu
	    all_n_spikes
	    
	    %PS results
	    PS_limit = min(all_n_spikes(2,:))
	    PS_tf_std = std(all_PS_tf(1:PS_limit,:),1,2);
	    
	    %RK results
	    RK_limit = min(all_n_spikes(1,:))
	    RK_tf_std = std(all_RK_tf(1:RK_limit,:),1,2);
	    
	    %BS results
	    BS_limit = min(all_n_spikes(3,:))
	    BS_tf_std = std(all_BS_tf(1:BS_limit,:),1,2);
	    
	%    keyboard
	%    pause
	%     
	%     continue
	    for cnd = 1:n_cnd
	      RK_self_ind = min(find(all_RK_nrn(:,cnd)-all_RK_nrn(:,end))-1);
	      RK_PS_ind = min(find(all_RK_nrn(:,cnd)-all_PS_nrn(:,end))-1);
	      RK_BS_ind = min(find(all_RK_nrn(:,cnd)-all_BS_nrn(:,end))-1);
	      if(isempty(RK_self_ind))
	        RK_self_acc(cnd) = t_end;
	      elseif(RK_self_ind == 0)
	        RK_self_acc(cnd) = 0;
	      else
	        RK_self_acc(cnd) = all_RK_tf(RK_self_ind,cnd);
	      end
	      if(isempty(RK_PS_ind))
	        RK_PS_acc(cnd) = t_end;
	      elseif(RK_PS_ind == 0)
	        RK_PS_acc(cnd) = 0;
	      else
	        RK_PS_acc(cnd) = all_RK_tf(RK_PS_ind,cnd);
	      end
	      if(isempty(RK_BS_ind))
	        RK_BS_acc(cnd) = t_end;
	      elseif(RK_BS_ind == 0)
	        RK_BS_acc(cnd) = 0;
	      else
	        RK_BS_acc(cnd) = all_RK_tf(RK_BS_ind,cnd);
	      end
	      
	      PS_self_ind = min(find(all_PS_nrn(:,cnd)-all_PS_nrn(:,end))-1);
	      PS_RK_ind = min(find(all_PS_nrn(:,cnd)-all_RK_nrn(:,end))-1);
	      PS_BS_ind = min(find(all_PS_nrn(:,cnd)-all_BS_nrn(:,end))-1);
	      if(isempty(PS_self_ind))
	        PS_self_acc(cnd) = t_end;
	      elseif(PS_self_ind == 0)
	        PS_self_acc(cnd) = 0;
	      else
	        PS_self_acc(cnd) = all_PS_tf(PS_self_ind,cnd);
	      end
	      if(isempty(PS_RK_ind))
	        PS_RK_acc(cnd) = t_end;
	      elseif(PS_RK_ind == 0)
	        PS_RK_acc(cnd) = 0;
	      else
	        PS_RK_acc(cnd) = all_PS_tf(PS_RK_ind,cnd);
	      end
	      if(isempty(PS_BS_ind))
	        PS_BS_acc(cnd) = t_end;
	      elseif(PS_BS_ind == 0)
	        PS_BS_acc(cnd) = 0;
	      else
	        PS_BS_acc(cnd) = all_PS_tf(PS_BS_ind,cnd);
	      end
	      
	      BS_self_ind = min(find(all_BS_nrn(:,cnd)-all_BS_nrn(:,end))-1);
	      BS_RK_ind = min(find(all_BS_nrn(:,cnd)-all_RK_nrn(:,end))-1);
	      BS_PS_ind = min(find(all_BS_nrn(:,cnd)-all_PS_nrn(:,end))-1);
	      if(isempty(BS_self_ind))
	        BS_self_acc(cnd) = t_end;
	      elseif(BS_self_ind == 0)
	        BS_self_acc(cnd) = 0;
	      else
	        BS_self_acc(cnd) = all_BS_tf(BS_self_ind,cnd);
	      end
	      if(isempty(BS_RK_ind))
	        BS_RK_acc(cnd) = t_end;
	      elseif(BS_RK_ind == 0)
	        BS_RK_acc(cnd) = 0;
	      else
	        BS_RK_acc(cnd) = all_BS_tf(BS_RK_ind,cnd);
	      end
	      if(isempty(BS_PS_ind))
	        BS_PS_acc(cnd) = t_end;
	      elseif(BS_PS_ind == 0)
	        BS_PS_acc(cnd) = 0;
	      else
	        BS_PS_acc(cnd) = all_BS_tf(BS_PS_ind,cnd);
	      end
	    end
	    acc = [RK_self_acc;RK_BS_acc;RK_PS_acc;BS_RK_acc;BS_self_acc;BS_PS_acc;PS_RK_acc;PS_BS_acc;PS_self_acc]
	    filename = [file_stub,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'.mat']
	    
	    save(filename,'ip','all_RK_tf','all_RK_nrn','all_PS_tf','all_PS_nrn','all_BS_tf','all_BS_nrn','all_t_cpu','all_n_spikes','all_RK_v','all_PS_v','all_BS_v','acc','ps_order_mean','ps_order_max','sub_mean','sub_max','bs_fails')
	  end
  end
end