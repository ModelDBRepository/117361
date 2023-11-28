% Analysis and display function for tm injection current simulations
% Written by Dr Robert Stewart for Stewart & Bair, 2009
function [] = tm_inj_analysis()
  disp('Injection current simulation analysis');
  %Set colours and marker sizes 
  rk_c = [1,0.6,0]; bs_c = [0,0.8,0]; ps_c = [1,0.0,1]; msize = 5;  
  
  load tm_inj_results_10_1.mat %Load one spike results
  t_RK1 = t_RK; t_BS1 = t_BS; t_PS1 = t_PS; %Store values
  nrm = sqrt(length(t_ref)); %normalisation for errorbars
  
  %Calculate maximum self-consistency value (comparing conditions 14 and 15)
  RK_self_err1 = mean(abs(all_RK_v(:,15)-all_RK_v(:,14)))
  BS_self_err1 = mean(abs(all_BS_v(:,15)-all_BS_v(:,14)))
  PS_self_err1 = mean(abs(all_PS_v(:,15)-all_PS_v(:,14)))
  
  RK_acc1 = 1./RK_err(1,:); BS_acc1 = 1./BS_err(1,:); PS_acc1 = 1./PS_err(1,:);
  RK_per1 = RK_acc1./mean(t_RK1);
	BS_per1 = BS_acc1./mean(t_BS1);
	PS_per1 = PS_acc1./mean(t_PS1);
  
  load tm_inj_results_10_10.mat %Load 10 spike results
  t_RK10 = t_RK; t_BS10 = t_BS; t_PS10 = t_PS; %Store values
  
  %Calculate maximum self-consistency value (comparing conditions 14 and 15)
  RK_self_err10 = mean(abs(all_RK_v(:,15)-all_RK_v(:,14)))
  BS_self_err10 = mean(abs(all_BS_v(:,15)-all_BS_v(:,14)))
  PS_self_err10 = mean(abs(all_PS_v(:,15)-all_PS_v(:,14)))
  
  RK_acc10 = 1./RK_err(1,:); BS_acc10 = 1./BS_err(1,:); PS_acc10 = 1./PS_err(1,:);
  RK_per10 = RK_acc10./mean(t_RK10);
	BS_per10 = BS_acc10./mean(t_BS10);
	PS_per10 = PS_acc10./mean(t_PS10);

  figure(1); clf; set(gcf,'Position',[400,300,350,450]); subplot(3,1,1)  
  semilogy(mean(t_RK1),'^-','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); 
  hold on  
  semilogy(mean(t_BS1),'o-','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
	semilogy(mean(t_PS1),'s-','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize); 
  semilogy(mean(t_RK10),'+:','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); 
  semilogy(mean(t_BS10),'x:','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
	semilogy(mean(t_PS10),'*:','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize);
  
  axis([0.7, 15.3 5*10^-2 10^4])
	set(gca,'XTick',1:15);
	set(gca,'XTickLabel',[]);
  set(gca,'Position',[0.17,0.72,.82,.26])
	
	set(gca,'YTick',[10^-1,10^0,10^1,10^2,10^3,10^4]);
	set(gca,'YTickLabel',{'0.1','1','10','100','1000','10000'});
	ylabel('Simulation time (s)')
  box on
  pos = get(get(gca,'ylabel'),'Position');
  pos(1) = pos(1)/1.5;
  pos(2) = pos(2)-4;
  set(get(gca,'ylabel'),'Position',pos)

  %Add custom legend
  semilogy(3.5,3500,'^','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([3:4],[3500,3500],'-','Color',rk_c); text(1,3500,'RK  1')
  
  semilogy(6.5,3500,'+','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([6:7],[3500,3500],':','Color',rk_c); text(5,3500,'10')

  semilogy(3.5,600,'o','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([3:4],[600,600],'-','Color',bs_c); text(1,600,'BS  1')
  
  semilogy(6.5,600,'x','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([6:7],[600,600],':','Color',bs_c); text(5,600,'10')
  
  semilogy(3.5,100,'s','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([3:4],[100,100],'-','Color',ps_c); text(1,100,'PS  1')
    
  semilogy(6.5,100,'*','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize);
  semilogy([6:7],[100,100],':','Color',ps_c); text(5,100,'10')
 
  subplot(3,1,2)
  semilogy(RK_acc1,'^-','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); 
  hold on  
  semilogy(BS_acc1,'o-','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
	semilogy(PS_acc1,'s-','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize); 
  semilogy(RK_acc10,'+:','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); 
  hold on  
  semilogy(BS_acc10,'x:','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
	semilogy(PS_acc10,'*:','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize); 
	set(gca,'XTick',1:15);
	set(gca,'XTickLabel',[]);
  set(gca,'Position',[0.17,0.41,.82,.26])
	axis([0.7, 15.3 10^0 10^15])
	set(gca,'YTick',[10^0,10^5,10^10,10^15]);
	ylabel('Accuracy (mV ^{-1})')
  box on
  
	subplot(3,1,3)
	semilogy(RK_per1,'^-','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); 
  hold on  
  semilogy(BS_per1,'o-','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
	semilogy(PS_per1,'s-','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize); 
  semilogy(RK_per10,'+:','Color',rk_c,'MarkerFaceColor','None','MarkerSize',msize); 
  hold on  
  semilogy(BS_per10,'+:','Color',bs_c,'MarkerFaceColor','None','MarkerSize',msize);
	semilogy(PS_per10,'*:','Color',ps_c,'MarkerFaceColor','None','MarkerSize',msize); 
	set(gca,'XTick',1:15);
  set(gca,'Position',[0.17,0.09,.82,.26])
	xlabel('Error tolerance condition')
	axis([0.7, 15.3 10^0 10^15])
  set(gca,'YTick',[10^0,10^5,10^10,10^15]);
  set(gca,'XTick',[1,5,10,15]);
	ylabel('Performance')
  box on
end