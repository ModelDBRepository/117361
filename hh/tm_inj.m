% Injection current simulation function for Traub-Miles neuron
% Written by Dr Robert Stewart for Stewart & Bair, 2009
function [] = tm_inj(n_spikes,tols,n_tols,dt_vals,ip,fp)
  if(n_spikes==1)
    I_inj = 20; fp(9) = I_inj; %20 pA -> 1 spike
  elseif(n_spikes==10)
    I_inj = 52; fp(9) = I_inj; %52 pA -> 10 spikes
  end

  for(in_seed=1:10) %Use in_seed as repeat index
    disp(['Input Condition: ',num2str(in_seed)]);
    ps_only = 0; %use 1 for debugging ps code only
    ip(4) = int32(in_seed); ip(5) = int32(ps_only);
    for(cnd = 1:n_tols)
      disp(['Tolerance Condition: ',num2str(cnd)]);
    	%Set error tolerance and RK dt given condition
      tol = tols(cnd); dt_rk = dt_vals(cnd);
      fp(10) = tol; fp(11) = dt_rk;
      [PS_v,RK_v,BS_v,t_cpu] = tm_ps(fp,ip);
      all_PS_v(:,cnd)=PS_v;
      all_BS_v(:,cnd)=BS_v;
      all_RK_v(:,cnd)=RK_v;
      t_PS(in_seed,cnd)=t_cpu(1);
      t_BS(in_seed,cnd)=t_cpu(2);
      t_RK(in_seed,cnd)=t_cpu(3);
    end
    cnd=n_tols+1;
    tol=0; ps_only=1; fp(10) = tol; ip(5) = ps_only;
    [PS_v,RK_v,BS_v,t_cpu] = tm_ps(fp,ip);
    ref_v=PS_v*ones(1,cnd-1);
    t_ref(in_seed) = t_cpu(3);
    
    RK_err(in_seed,:) = mean(abs(all_RK_v-ref_v));
    BS_err(in_seed,:) = mean(abs(all_BS_v-ref_v));
    PS_err(in_seed,:) = mean(abs(all_PS_v-ref_v));  
    
    RK_per(in_seed,:) = 1./(RK_err(in_seed,:).*t_RK(in_seed,:))
    BS_per(in_seed,:) = 1./(BS_err(in_seed,:).*t_BS(in_seed,:))
    PS_per(in_seed,:) = 1./(PS_err(in_seed,:).*t_PS(in_seed,:))
  end
  n_neurons = ip(6);
  filename = ['tm_inj_results_',num2str(n_neurons),'_',num2str(n_spikes)]
  save(filename);
end