% Analysis and display function for recurrent network simulations
% Written by Dr Robert Stewart for Stewart & Bair, 2009
function [] = iz_fix_analysis(fp)
  disp('Running analysis...');
  n_neurons = 4000; t_end = 1000; trace_ind = 5; %4 in c code -> 5 here
  max_spikes = 200000; vr = fp(6);
  syn_seed = 5; in_seed = 9; %Vary these to view results from different sims
  file_stub = 'iz_fix';
  filename = [file_stub,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'.mat']
	load(filename);
	
	%Get objective divergence points to annotate graph
  dval = 1;
	ref_vals = all_PS_v(:,4)*ones(1,3);
	rk_dif = abs(all_RK_v(:,1:3)-ref_vals);
	bs_dif = abs(all_BS_v(:,1:3)-ref_vals);
	ps_dif = abs(all_PS_v(:,1:3)-ref_vals);
	div_rk(1) = min(find(rk_dif(:,1)>dval))-1;
  div_rk(2) = min(find(rk_dif(:,2)>dval))-1;
  div_rk(3) = min(find(rk_dif(:,3)>dval))-1;
  
  div_bs(1) = min(find(bs_dif(:,1)>dval))-1;
  div_bs(2) = min(find(bs_dif(:,2)>dval))-1;
  div_bs(3) = min(find(bs_dif(:,3)>dval))-1;
  
  div_ps(1) = min(find(ps_dif(:,1)>dval))-1;
  div_ps(2) = min(find(ps_dif(:,2)>dval))-1;
  temp = min(find(ps_dif(:,3)>dval))-1;
  if(~isempty(temp))
  	div_ps(3) = temp;
  else
  	disp('Perfect PS');
  end
  
  %Get reference solution
  ref_cnd = 4; file_stub2 = 'iz_bench_fix';
  ref_filename = [file_stub2,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'_',num2str(ref_cnd-1),'_ps_trace']
  fid_ref = fopen(ref_filename); ref_trace = fread(fid_ref,[6,inf],'double'); fclose(fid_ref);
  
  %Set colours
  rk_c = [1,0.6,0]; bs_c = [0,0.8,0]; ps_c = [1,0.0,1];  
  
  %Trace figure
  figure(1); clf;
  set(gcf,'Position',[800,600,350,400])
  for cnd = 1:3
  	subplot(3,1,cnd)
  	filename = [file_stub2,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'_',num2str(cnd-1),'_rk_trace'];
  	fid = fopen(filename); rk_trace = fread(fid,[5,inf],'double'); fclose(fid);  	
  	plot(rk_trace(1,:),rk_trace(2,:)+vr,'Color',rk_c); hold on
  	filename = [file_stub2,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'_',num2str(cnd-1),'_bs_trace'];
  	fid = fopen(filename); bs_trace = fread(fid,[6,inf],'double'); fclose(fid);  	
  	plot(bs_trace(1,:),bs_trace(2,:)+vr,'Color',bs_c);  	
  	filename = [file_stub2,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'_',num2str(cnd-1),'_ps_trace'];
  	fid = fopen(filename); ps_trace = fread(fid,[6,inf],'double'); fclose(fid);
  	plot(ps_trace(1,:),ps_trace(2,:)+vr,'Color',ps_c); 
  	plot(ref_trace(1,:),ref_trace(2,:)+vr,'k-')
  	axis([0 t_end -80 0])
  	set(gca,'YTick',[-80,-40,0])
  	set(gca,'Position',[0.12, 0.09+(3-cnd)*.31, .84, .26])

  	if(cnd==1)
  		legend('RK','BS','PS','ref','Location','East')
  	end
  	if(cnd==2)
  		ylabel('Membrane potential (mV)'); 
  		%set(get(gca,'ylabel'),'Position',[-83.3333 0 1.00011])
  	end
  	if(cnd<3)
  		set(gca,'XTickLabel',[])
  	else
  		xlabel('Time (ms)');
  	end
  	
  	%Annotate with objective divergence points  	
  	plot([div_rk(cnd),div_rk(cnd)],[-80 -1],'^-','Color',rk_c,'MarkerFaceColor',rk_c,'MarkerSize',8)
	  plot([div_bs(cnd),div_bs(cnd)],[-80,-1],'o-','Color',bs_c,'MarkerFaceColor',bs_c,'MarkerSize',8)
	  if(cnd<3)
	  	plot([div_ps(cnd),div_ps(cnd)],[-80,-1],'s-','Color',ps_c,'MarkerFaceColor',ps_c,'MarkerSize',8)
	  end
  end
  %Annotate with condition numbers
  annotation('textbox',[0.97,0.93,0.03,0.03],'EdgeColor',[1 1 1],'String','1','FontWeight','bold','Margin', 0,'FitBoxToText','on');
  annotation('textbox',[0.97,0.63,0.03,0.03],'EdgeColor',[1 1 1],'String','2','FontWeight','bold','Margin', 0,'FitBoxToText','on');
  annotation('textbox',[0.97,0.32,0.03,0.03],'EdgeColor',[1 1 1],'String','3','FontWeight','bold','Margin', 0,'FitBoxToText','on');
 
 	%Population raster plot
 	figure(2); clf; 
 	dotsize = 12;
 	set(gcf,'Position',[1200,600,350,400])
 	for cnd = 1:3
 		div_rk_raster = 1000;
  	div_bs_raster = 1000;
  	div_ps_raster = 1000;
  	subplot(3,1,cnd); hold on
  	box on
  	set(gca,'Position',[0.12, 0.09+(3-cnd)*.31, .84, .26])
  	if(cnd<3)
  		set(gca,'XTickLabel',[])
  	else
  		xlabel('Time (ms)');
  	end
  	axis([0 1000 0 15.9])
  	for nrn_ind = 15:-1:1 %reverse order for legend
	    ind = find(all_RK_nrn(:,cnd)==nrn_ind); t1 = all_RK_tf(ind,cnd);
	    if(~isempty(ind))
	      plot(t1,ones(size(t1))*nrn_ind,'.','Color',rk_c,'Markersize',dotsize);
      end	    
	    ind = find(all_BS_nrn(:,cnd)==nrn_ind); t2 = all_BS_tf(ind,cnd);
	    if(~isempty(ind))
	      plot(t2,ones(size(t2))*nrn_ind,'.','Color',bs_c,'Markersize',dotsize);
      end	    
	    ind = find(all_PS_nrn(:,cnd)==nrn_ind); t3 = all_PS_tf(ind,cnd);
	    if(~isempty(ind))
	      plot(t3,ones(size(t3))*nrn_ind,'.','Color',ps_c,'Markersize',dotsize);
      end	    
	    ind = find(all_PS_nrn(:,4)==nrn_ind); t4 = all_PS_tf(ind,4);
	    if(~isempty(ind))
	      plot(t4,ones(size(t4))*nrn_ind,'k.','Markersize',dotsize);
	    end
	    
	    %Get divergence and annotate
	    div = 1000;
	    for i = 1:length(t4) %loop over reference spikes
	    	if(i>length(t1))
	    		div = t4(i); break
	    	end
	    	if(abs(t4(i)-t1(i))>dval)
	    		div = min([t4(i),t1(i)]); break
	    	end
	    end
	    if(div<div_rk_raster)
	    	div_rk_raster=div;
	    end
	    
	    div = 1000;
	    for i = 1:length(t4) %loop over reference spikes
	    	if(i>length(t2))
	    		div = t4(i); break
	    	end
	    	if(abs(t4(i)-t2(i))>dval)
	    		div = min([t4(i),t2(i)]); break
	    	end
	    end
	    if(div<div_bs_raster) 
	    	div_bs_raster=div;
	    end
	    
	    div = 1000;
	    for i = 1:length(t4) %loop over reference spikes
	    	if(i>length(t3))
	    		div = t4(i); break
	    	end
	    	if(abs(t4(i)-t3(i))>dval)
	    		div = min([t4(i),t3(i)]); break
	    	end
	    end
	    if(div<div_ps_raster) 
	    	div_ps_raster=div;
	    end
	    if(nrn_ind==15)
	  		if(cnd==1)
  				legend('RK','BS','PS','ref','Location','East')
  			end
  		end
    end	  
	  if(cnd==2)
  		ylabel('Neuron index');
  	end
	  div_rk_rasters(cnd) = div_rk_raster;
	  div_bs_rasters(cnd) = div_bs_raster;
	  div_ps_rasters(cnd) = div_ps_raster;
	  plot([div_rk_raster,div_rk_raster],[0,15.8],'^-','Color',rk_c,'MarkerFaceColor',rk_c,'MarkerSize',8)
	  plot([div_bs_raster,div_bs_raster],[0,15.8],'o-','Color',bs_c,'MarkerFaceColor',bs_c,'MarkerSize',8)
	  if(cnd<3)
	  	plot([div_ps_raster,div_ps_raster],[0,15.8],'s-','Color',ps_c,'MarkerFaceColor',ps_c,'MarkerSize',8)
	  end
  end
  %Annotate with condition numbers  
  annotation('textbox',[0.97,0.93,0.03,0.03],'EdgeColor',[1 1 1],'String','1','FontWeight','bold','Margin', 0,'FitBoxToText','on');
  annotation('textbox',[0.97,0.63,0.03,0.03],'EdgeColor',[1 1 1],'String','2','FontWeight','bold','Margin', 0,'FitBoxToText','on');
  annotation('textbox',[0.97,0.32,0.03,0.03],'EdgeColor',[1 1 1],'String','3','FontWeight','bold','Margin', 0,'FitBoxToText','on');
  
  disp('Gathering time and divergence data');
  syn_seeds = 1:5; in_seeds = 1:10; ind = 0;
  for syn_seed=syn_seeds
  	for in_seed = in_seeds
  		filename = [file_stub,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'.mat']
			load(filename);	ind=ind+1;
			%Gather CPU Times
			RK_t_cpu(ind,:)=all_t_cpu(1,:);BS_t_cpu(ind,:)=all_t_cpu(2,:);PS_t_cpu(ind,:)=all_t_cpu(3,:);
      total_n_spikes(ind,:) = all_n_spikes(all_n_spikes>0);
			
			%adaptive processing and failures
  		all_bs_fails(ind,:) = bs_fails;
  		all_ps_order_mean(ind,:) = ps_order_mean;	all_ps_order_max(ind,:) = ps_order_max;
  		all_sub_mean(ind,:) = sub_mean;	all_sub_max(ind,:) = sub_max;	
  		
			%Get reference spike sequence
			ref = all_PS_nrn(:,4);
			for cnd = 1:3
				RK_ind = min(find(all_RK_nrn(:,cnd)-ref)-1);
      	if(isempty(RK_ind))
        	all_RK_acc(ind,cnd) = t_end;
      	elseif(RK_ind == 0)
        	all_RK_acc(ind,cnd) = 0;
      	else
        	all_RK_acc(ind,cnd) = all_RK_tf(RK_ind,cnd);
      	end
      end			
			for cnd = 1:3
				BS_ind = min(find(all_BS_nrn(:,cnd)-ref)-1);
      	if(isempty(BS_ind))
        	all_BS_acc(ind,cnd) = t_end;
      	elseif(BS_ind == 0)
        	all_BS_acc(ind,cnd) = 0;
      	else
        	all_BS_acc(ind,cnd) = all_BS_tf(BS_ind,cnd);
      	end
      end			
			for cnd = 1:3
				PS_ind = min(find(all_PS_nrn(:,cnd)-ref)-1);
      	if(isempty(PS_ind))
        	all_PS_acc(ind,cnd) = t_end;
      	elseif(PS_ind == 0)
        	all_PS_acc(ind,cnd) = 0;
      	else
        	all_PS_acc(ind,cnd) = all_PS_tf(PS_ind,cnd);
      	end
			end
		end
  end
  
  %Figure 4 for paper putting together time, accuracy, performance
  figure(4); clf;
	nrm = sqrt(length(RK_t_cpu));
	set(gcf,'Position',[800,300,350,450])
	%TIME
	subplot(3,1,1)
 	semilogy(mean(RK_t_cpu(:,1:3)),'^-','Color',rk_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5); hold on	
	semilogy(mean(BS_t_cpu(:,1:3)),'o--','Color',bs_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5);
	semilogy(mean(PS_t_cpu(:,1:3)),'s:','Color',ps_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5);
	legend('RK','BS','PS','ref','Location','NorthWest')
	axis([0.9, 3.1 10^0 10^4])
	set(gca,'XTick',[1,2,3]);
	set(gca,'XTickLabel',[]);
	ylabel('Simulation time (s)');  
  set(gca,'YTick',[10^0,10^1,10^2,10^3,10^4]);	
	set(gca,'YTickLabel',[1,10,100,1000,10000]);
	set(gca,'Position',[0.15,0.72,.84,.26])
  p = get(get(gca,'ylabel'),'Position')
  p(1) = p(1)+.05;
  set(get(gca,'ylabel'),'Position',p)
  
  all_RK_acc = all_RK_acc/1000;all_BS_acc = all_BS_acc/1000;all_PS_acc = all_PS_acc/1000;
  
	%ACCURACY
	subplot(3,1,2)
 	plot(mean(all_RK_acc),'^-','Color',rk_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5);	hold on	
	plot(mean(all_BS_acc),'o--','Color',bs_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5);
	plot(mean(all_PS_acc),'s:','Color',ps_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5);
	set(gca,'XTick',[1,2,3,]);
	set(gca,'XTickLabel',[]);
	ylabel({'Duration of';'agreement (s)'})
	axis([0.9, 3.1 0 1])  
	set(gca,'Position',[0.15,0.41,.84,.26])
  
  %PERFORMANCE (ACC/TIME)
	RK_perf = all_RK_acc./RK_t_cpu(:,1:3);
	BS_perf = all_BS_acc./BS_t_cpu(:,1:3);
	PS_perf = all_PS_acc./PS_t_cpu(:,1:3);
	subplot(3,1,3)
 	plot(mean(RK_perf),'^-','Color',rk_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5); hold on	
	plot(mean(BS_perf),'o--','Color',bs_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5);
	plot(mean(PS_perf),'s:','Color',ps_c,'MarkerFaceColor','None','MarkerSize',8,'LineWidth',1.5);
	xlabel('Error tolerance condition')
	set(gca,'XTick',[1,2,3]);
	axis([0.9, 3.1 0 .05])
	ylabel('Performance');  
	set(gca,'Position',[0.15,0.07,.84,.28])
  
  %Gather simulation time data from smaller simulations for comparison
  disp('Gathering time and divergence data (1000 cells)');
  n_neurons = 1000; syn_seeds = 1:5; in_seeds = 1:10; ind = 0;
  for syn_seed=syn_seeds
  	for in_seed = in_seeds
  		filename = [file_stub,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'.mat']
			load(filename);	ind=ind+1;
			%Gather CPU Times
			RK_t_cpu_1000(ind,:)=all_t_cpu(1,:);BS_t_cpu_1000(ind,:)=all_t_cpu(2,:);PS_t_cpu_1000(ind,:)=all_t_cpu(3,:);
			total_n_spikes_1000(ind,:) = all_n_spikes(all_n_spikes>0);
			
			%adaptive processing and failures
  		all_bs_fails_1000(ind,:) = bs_fails;
  		all_ps_order_mean_1000(ind,:) = ps_order_mean; all_ps_order_max_1000(ind,:) = ps_order_max;
  		all_sub_mean_1000(ind,:) = sub_mean; all_sub_max_1000(ind,:) = sub_max;	
  		
			ref = all_PS_nrn(:,4); %Get reference spike sequence
			for cnd = 1:3
				RK_ind = min(find(all_RK_nrn(:,cnd)-ref)-1);
      	if(isempty(RK_ind))
        	all_RK_acc_1000(ind,cnd) = t_end;
      	elseif(RK_ind == 0)
        	all_RK_acc_1000(ind,cnd) = 0;
      	else
        	all_RK_acc_1000(ind,cnd) = all_RK_tf(RK_ind,cnd);
      	end
      end			
			for cnd = 1:3
				BS_ind = min(find(all_BS_nrn(:,cnd)-ref)-1);
      	if(isempty(BS_ind))
        	all_BS_acc_1000(ind,cnd) = t_end;
      	elseif(BS_ind == 0)
        	all_BS_acc_1000(ind,cnd) = 0;
      	else
        	all_BS_acc_1000(ind,cnd) = all_BS_tf(BS_ind,cnd);
      	end
      end			
			for cnd = 1:3
				PS_ind = min(find(all_PS_nrn(:,cnd)-ref)-1);
      	if(isempty(PS_ind))
        	all_PS_acc_1000(ind,cnd) = t_end;
      	elseif(PS_ind == 0)
        	all_PS_acc_1000(ind,cnd) = 0;
      	else
        	all_PS_acc_1000(ind,cnd) = all_PS_tf(PS_ind,cnd);
      	end
			end
		end
  end
  
  disp('Gathering time and divergence data (2000 cells)');
  n_neurons = 2000; syn_seeds = 1:5; in_seeds = 1:10; ind = 0;
  for syn_seed=syn_seeds
  	%disp(['Config seed = ',num2str(syn_seed)])
  	for in_seed = in_seeds
    	%disp(['Input seed = ',num2str(in_seed)])
  		filename = [file_stub,'_',num2str(n_neurons),'_',num2str(syn_seed),'_',num2str(in_seed),'.mat']
			load(filename);	ind=ind+1;
			%Gather CPU Times
			RK_t_cpu_2000(ind,:)=all_t_cpu(1,:);BS_t_cpu_2000(ind,:)=all_t_cpu(2,:);PS_t_cpu_2000(ind,:)=all_t_cpu(3,:);
			total_n_spikes_2000(ind,:) = all_n_spikes(all_n_spikes>0);
			
			%adaptive processing and failures
  		all_bs_fails_2000(ind,:) = bs_fails;
  		all_ps_order_mean_2000(ind,:) = ps_order_mean; all_ps_order_max_2000(ind,:) = ps_order_max;
  		all_sub_mean_2000(ind,:) = sub_mean; all_sub_max_2000(ind,:) = sub_max;	
  		
			%Get reference spike sequence
			ref = all_PS_nrn(:,4);
			for cnd = 1:3
				RK_ind = min(find(all_RK_nrn(:,cnd)-ref)-1);
      	if(isempty(RK_ind))
        	all_RK_acc_2000(ind,cnd) = t_end;
      	elseif(RK_ind == 0)
        	all_RK_acc_2000(ind,cnd) = 0;
      	else
        	all_RK_acc_2000(ind,cnd) = all_RK_tf(RK_ind,cnd);
      	end
      end			
			for cnd = 1:3
				BS_ind = min(find(all_BS_nrn(:,cnd)-ref)-1);
      	if(isempty(BS_ind))
        	all_BS_acc_2000(ind,cnd) = t_end;
      	elseif(BS_ind == 0)
        	all_BS_acc_2000(ind,cnd) = 0;
      	else
        	all_BS_acc_2000(ind,cnd) = all_BS_tf(BS_ind,cnd);
      	end
      end			
			for cnd = 1:3
				PS_ind = min(find(all_PS_nrn(:,cnd)-ref)-1);
      	if(isempty(PS_ind))
        	all_PS_acc_2000(ind,cnd) = t_end;
      	elseif(PS_ind == 0)
        	all_PS_acc_2000(ind,cnd) = 0;
      	else
        	all_PS_acc_2000(ind,cnd) = all_PS_tf(PS_ind,cnd);
      	end
			end
		end
  end
  %Get time ratios
  mean_BS_t_ratios(1) = mean(mean(BS_t_cpu_2000(:,1:3))./mean(BS_t_cpu_1000(:,1:3)));
  std_BS_t_ratios(1) = std(mean(BS_t_cpu_2000(:,1:3))./mean(BS_t_cpu_1000(:,1:3)));
  mean_BS_t_ratios(2) = mean(mean(BS_t_cpu(:,1:3))./mean(BS_t_cpu_2000(:,1:3)))
  std_BS_t_ratios(2) = std(mean(BS_t_cpu(:,1:3))./mean(BS_t_cpu_2000(:,1:3)))
  
  mean_PS_t_ratios(1) = mean(mean(PS_t_cpu_2000(:,1:3))./mean(PS_t_cpu_1000(:,1:3)));
  std_PS_t_ratios(1) = std(mean(PS_t_cpu_2000(:,1:3))./mean(PS_t_cpu_1000(:,1:3)));
  mean_PS_t_ratios(2) = mean(mean(PS_t_cpu(:,1:3))./mean(PS_t_cpu_2000(:,1:3)))
  std_PS_t_ratios(2) = std(mean(PS_t_cpu(:,1:3))./mean(PS_t_cpu_2000(:,1:3)))
  
  mean_RK_t_ratios(1) = mean(mean(RK_t_cpu_2000(:,1:3))./mean(RK_t_cpu_1000(:,1:3)));
  std_RK_t_ratios(1) = std(mean(RK_t_cpu_2000(:,1:3))./mean(RK_t_cpu_1000(:,1:3)));
  mean_RK_t_ratios(2) = mean(mean(RK_t_cpu(:,1:3))./mean(RK_t_cpu_2000(:,1:3)))
  std_RK_t_ratios(2) = std(mean(RK_t_cpu(:,1:3))./mean(RK_t_cpu_2000(:,1:3)))

  load inj_results_final.mat %Get injection current results to compare
	temp = mean(t_BS10); t_BS_div = temp([1,9,15])
	temp = mean(t_PS10); t_PS_div = temp([1,9,15])
	temp = mean(t_RK10); t_RK_div = temp([1,9,15])
	
	BS_t_ratios_basic = mean(BS_t_cpu(:,1:3))./t_BS_div
	PS_t_ratios_basic = mean(PS_t_cpu(:,1:3))./t_PS_div
	RK_t_ratios_basic = mean(RK_t_cpu(:,1:3))./t_RK_div
	  
	BS_t_ratios_small = mean(BS_t_cpu_1000(1:10,1:3))./t_BS_div
	PS_t_ratios_small = mean(PS_t_cpu_1000(1:10,1:3))./t_PS_div
	RK_t_ratios_small = mean(RK_t_cpu_1000(1:10,1:3))./t_RK_div
end