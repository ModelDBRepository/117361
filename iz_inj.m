% Injection current simulation function
% Written by Dr Robert Stewart for Stewart & Bair, 2009
function [] = iz_inj(n_spikes,tols,n_tols,dt_vals,ip,fp)
  if(n_spikes==1)
    I_inj = 21; fp(1) = I_inj; %21 pA -> 1 spike
  elseif(n_spikes==10)
    I_inj = 30; fp(1) = I_inj; %30 pA -> 10 spikes
  end
  for(in_seed=1:10) %Use in_seed as repeat index
    in_seed
    ps_only = 0;
    ip(5) = int32(in_seed); ip(7) = int32(ps_only);
    for(cnd = 1:n_tols)	
      disp(['Condition: ',num2str(cnd)])
      ip(6) = int32(cnd-1);
      %Set error tolerance and RK dt given condition
      tol = tols(cnd); dt_rk = dt_vals(cnd); fp(2:3) = [tol,dt_rk];      
      [RK_tf,RK_nrn,PS_tf,PS_nrn,BS_tf,BS_nrn,t_cpu,RK_v,PS_v,BS_v,i_stats,f_stats] = iz_ps(fp,ip);
      RK_spike_times = unique(RK_tf(find(RK_tf)));
      BS_spike_times = unique(BS_tf(find(BS_tf)));
      PS_spike_times = unique(PS_tf(find(PS_tf)));
      all_spike_times = [RK_spike_times,BS_spike_times,PS_spike_times];
      all_PS_v(:,cnd)=PS_v; all_BS_v(:,cnd)=BS_v; all_RK_v(:,cnd)=RK_v;
      
      t_RK(in_seed,cnd)=t_cpu(1);t_BS(in_seed,cnd)=t_cpu(2);t_PS(in_seed,cnd)=t_cpu(3);  
      
      PS_mean_order(in_seed,cnd) = f_stats(1);
      PS_max_order(in_seed,cnd) = i_stats(1);
      BS_mean_crossings(in_seed,cnd) = f_stats(2)/2;
      BS_max_crossings(in_seed,cnd) = i_stats(2)/2;
      BS_fails(in_seed,cnd) = i_stats(3);
    end
    
    %run reference PS condition with zero tolerance
    cnd=n_tols+1; ip(6) = int32(cnd-1); ps_only=1; ip(7) = int32(ps_only); fp(2) = 0;
    [RK_tf,RK_nrn,PS_tf,PS_nrn,BS_tf,BS_nrn,t_cpu,RK_v,PS_v,BS_v,i_stats,f_stats] = iz_ps(fp,ip);
    ref_v=PS_v*ones(1,cnd-1);
    
    RK_err(in_seed,:) = mean(abs(all_RK_v-ref_v));
    BS_err(in_seed,:) = mean(abs(all_BS_v-ref_v));
    PS_err(in_seed,:) = mean(abs(all_PS_v-ref_v));  
    
    RK_per(in_seed,:) = 1./(RK_err(in_seed,:).*t_RK(in_seed,:))
    BS_per(in_seed,:) = 1./(BS_err(in_seed,:).*t_BS(in_seed,:))
    PS_per(in_seed,:) = 1./(PS_err(in_seed,:).*t_PS(in_seed,:))
    
    t_ref(in_seed) = t_cpu(3);
    ref_mean_order(in_seed) = f_stats(1); ref_max_order(in_seed) = i_stats(1);
  end
  n_neurons = ip(1);
  filename = ['inj_results_',num2str(n_neurons),'_',num2str(n_spikes)]
  save(filename);
end