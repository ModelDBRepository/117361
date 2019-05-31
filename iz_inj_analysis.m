% Analysis and display function for injection current simulations
% Written by Dr Robert Stewart for Stewart & Bair, 2009
function [] = iz_inj_analysis()
  disp('Injection current simulation analysis');
  %Set colours and marker sizes 
  rk_c = [1,0.6,0]; bs_c = [0,0.8,0]; ps_c = [1,0.0,1]; msize = 5;
  
  load inj_results_1000_1.mat %One Spike results
  mean_t_ref1 = mean(t_ref); std_t_ref1 = std(t_ref);
  nrm = sqrt(length(t_ref));  %normalisation for errorbars
   
  %Store values for one spike results
  t_RK1 = t_RK; t_BS1 = t_BS; t_PS1 = t_PS;
  BS_mean_crossings1 = BS_mean_crossings;
  PS_mean_order1 = PS_mean_order; PS_max_order1 = PS_max_order;
  
  RK_test1 = mean(t_RK1).*dt_vals; %time cost per step in ms
  PS_test1 = mean(t_PS1)./PS_mean_order1(1,:).*dt_vals(1); %time/step/order, ms
	BS_test1 = mean(t_BS1)./BS_mean_crossings1(1,:).*dt_vals(1);%time/step/crossing, ms
  
  %Self-consistency test (comparing conditions 14 and 15 only)
  RK_self_err1 = mean(abs(all_RK_v(:,15)-all_RK_v(:,14)))
  BS_self_err1 = mean(abs(all_BS_v(:,15)-all_BS_v(:,14)))
  PS_self_err1 = mean(abs(all_PS_v(:,15)-all_PS_v(:,14)))
  
  RK_acc1 = 1./RK_err(1,:); BS_acc1 = 1./BS_err(1,:); PS_acc1 = 1./PS_err(1,:);
  RK_per1 = RK_acc1./mean(t_RK1);
	BS_per1 = BS_acc1./mean(t_BS1);
	PS_per1 = PS_acc1./mean(t_PS1);
  
  load inj_results_1000_10.mat %10 Spike results
  mean_t_ref10 = mean(t_ref); std_t_ref10 = std(t_ref);
  
  %Store values for ten spike results
  t_RK10 = t_RK; t_BS10 = t_BS; t_PS10 = t_PS;
  BS_mean_crossings10 = BS_mean_crossings;
  PS_mean_order10 = PS_mean_order; PS_max_order10 = PS_max_order;
  
  RK_test10 = mean(t_RK10).*dt_vals; %time cost per step in ms
  PS_test10 = mean(t_PS10)./PS_mean_order10(1,:).*dt_vals(1); %time/step/order
	BS_test10 = mean(t_BS10)./BS_mean_crossings10(1,:).*dt_vals(1); %time/step/crossing
  
  %Self-consistency test (comparing conditions 14 and 15 only)
  RK_self_err10 = mean(abs(all_RK_v(:,15)-all_RK_v(:,14)))
  BS_self_err10 = mean(abs(all_BS_v(:,15)-all_BS_v(:,14)))
  PS_self_err10 = mean(abs(all_PS_v(:,15)-all_PS_v(:,14)))
  
  RK_acc10 = 1./RK_err(1,:); BS_acc10 = 1./BS_err(1,:); PS_acc10 = 1./PS_err(1,:);
  RK_per10 = RK_acc10./mean(t_RK10);
	BS_per10 = BS_acc10./mean(t_BS10);
	PS_per10 = PS_acc10./mean(t_PS10);
  
  save inj_results_final.mat %save final results to file

  figure(1); clf;
  set(gcf,'Position',[400,300,350,450])
  subplot(3,1,1)
  semilogy(mean(t_RK1),'^-','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); hold on  
  semilogy(mean(t_BS1),'o-','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy(mean(t_PS1),'s-','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize); 
  semilogy(mean(t_RK10),'+:','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); 
  semilogy(mean(t_BS10),'x:','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy(mean(t_PS10),'*:','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize);
  axis([0.7, 15.3 2^-1 10^3])

  %Manually add custom legend
  semilogy(3.5,350,'^','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([3:4],[350,350],'-','Color',rk_c); text(1,350,'RK  1')

  semilogy(6.5,350,'+','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([6:7],[350,350],':','Color',rk_c); text(5,350,'10')

  semilogy(3.5,100,'o','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([3:4],[100,100],'-','Color',bs_c); text(1,100,'BS  1')

  semilogy(6.5,100,'x','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([6:7],[100,100],':','Color',bs_c); text(5,100,'10')

  semilogy(3.5,30,'s','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([3:4],[30,30],'-','Color',ps_c); text(1,30,'PS  1')

  semilogy(6.5,30,'*','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([6:7],[30,30],':','Color',ps_c); text(5,30,'10')

  set(gca,'XTick',1:15);
  set(gca,'XTickLabel',[]);
  set(gca,'Position',[0.15,0.72,.84,.26])

  set(gca,'YTick',[10^0,10^1,10^2,10^3]);
  set(gca,'YTickLabel',{'1','10','100','1000'});
  ylabel('Simulation time (s)')
  box on
  pos = get(get(gca,'ylabel'),'Position');
  pos(1) = pos(1)/1.5;
  pos(2) = pos(2)-4;
  set(get(gca,'ylabel'),'Position',pos)

  %adaptive processing stats
  subplot(3,1,2) 
  plot(BS_mean_crossings1(1,:),'o-','Color',bs_c,'MarkerSize',msize); hold on
  plot(PS_mean_order1(1,:),'s-','Color',ps_c,'MarkerSize',msize);
  plot(PS_max_order1(1,:),'k-d','MarkerSize',msize);

  plot(BS_mean_crossings10(1,:),'x:','Color',bs_c,'MarkerSize',msize); hold on
  plot(PS_mean_order10(1,:),'*:','Color',ps_c,'MarkerSize',msize);
  plot(PS_max_order10(1,:),'k+:','MarkerSize',msize);

  text(1,43,'Mean BS crossings')
  text(1,35,'Mean PS order')
  text(1,27,'Max PS order')

  text(7,43,'1')
  text(7,35,'1')
  text(7,27,'1')

  plot(8,43,'o','Color',bs_c,'MarkerSize',msize);
  plot(7.5:8.5,[43,43],'-','Color',bs_c);

  plot(8,35,'s','Color',ps_c,'MarkerSize',msize);
  plot(7.5:8.5,[35,35],'-','Color',ps_c);

  plot(8,27,'kd','MarkerSize',msize);
  plot(7.5:8.5,[27,27],'k-');

  text(9.5,43,'10')
  text(9.5,35,'10')
  text(9.5,27,'10')

  plot(11,43,'x','Color',bs_c,'MarkerSize',msize);
  plot(10.5:11.5,[43,43],':','Color',bs_c);

  plot(11,35,'*','Color',ps_c,'MarkerSize',msize);
  plot(10.5:11.5,[35,35],':','Color',ps_c);

  plot(11,27,'k+','MarkerSize',msize);
  plot(10.5:11.5,[27,27],'k:');

  ylabel({'PS order /';'BS crossings'})
  axis([0.7, 15.3 0 50])
  set(gca,'YTick',[0,10,20,30 40 50]);
  set(gca,'XTick',1:15);
  set(gca,'XTickLabel',[]);
  set(gca,'Position',[0.15,0.41,.84,.26])

  subplot(3,1,3)
  semilogy(RK_acc1,'^-','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); hold on  
  semilogy(BS_acc1,'o-','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy(PS_acc1,'s-','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize); 
  semilogy(RK_acc10,'+:','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); hold on  
  semilogy(BS_acc10,'x:','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy(PS_acc10,'*:','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize); 
  set(gca,'Position',[0.15,0.07,.84,.28])
  axis([0.7, 15.3 10^0 10^20])
  set(gca,'YTick',[10^0,10^5,10^10,10^15,10^20]);
  ylabel('Accuracy (mV ^{-1})')
  xlabel('Error tolerance condition')
  set(gca,'YTick',[10^0,10^5,10^10,10^15,10^20]);
  set(gca,'XTick',[1,5,10,15]);
  box on
end