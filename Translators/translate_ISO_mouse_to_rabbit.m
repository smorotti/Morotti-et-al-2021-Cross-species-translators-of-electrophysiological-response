%% Translate the effects of ISO administration from mouse to rabbit (1 Hz)

clear
close all
clc

%% Select number of features to be used in mouse (2 or 4)
n_features_mouse = 4; % 2 or 4

if n_features_mouse == 2
    translator_MtoR_2 % APD90, APD50
elseif n_features_mouse == 4
    translator_MtoR_4 % APD90, APD50, CaTamp, CaTtau
end

%% Plotting options
color_human = [0 0.45 0.74];
color_rabbit = [0.85 0.33 0.1];
color_mouse = [0.47 0.67 0.19];

%% Rabbit - Hegyi et al. 2021
rabbit_control = [221.5581439 170.1345514 10 10 10 10]; % control rabbit 1 Hz
rabbit_control_sem = [8.001801947 7.371768043 1 1 1 1];
    rabbit_control_std = rabbit_control_sem*sqrt(18);

rabbit_iso = [188.19817 150.5785142 10 10 10 10]; % iso rabbit 1 Hz
rabbit_iso_sem = [5.262625723 6.6650293 1 1 1 1];
    rabbit_iso_std = rabbit_iso_sem*sqrt(6);

% *******************************************************************
figure,set(gcf,'color','w')
subplot(2,3,1),hold on
bar([rabbit_control(1);rabbit_iso(1)]','FaceColor',color_rabbit)
errorbar([1; 2]',[rabbit_control(1);rabbit_iso(1)]',[rabbit_control_std(1);rabbit_iso_std(1)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD90 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});

subplot(2,3,2),hold on
bar([rabbit_control(2);rabbit_iso(2)]','FaceColor',color_rabbit)
errorbar([1; 2]',[rabbit_control(2);rabbit_iso(2)]',[rabbit_control_std(2);rabbit_iso_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD50 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});
title('Rabbit (1-Hz) from Hegyi et al. 2021')
% *******************************************************************

%% Mouse - Edwards et al. 2014

% *******************************************************************
% WT - n = 17;
APD50 = 3.689505882; APD50_std = 1.890940989; APD50_sem = 0.458620555;
APD90 = 31.66589412; APD90_std = 15.4112023; APD90_sem = 3.737765582;
delta_Ca = 0.2054; delta_Ca_std = 0.133654606; delta_Ca_sem = 0.032416003;
Ca_tau = 259.9954471; Ca_tau_std = 106.1922159; Ca_tau_sem = 25.75539545;

% WT + ISO - n = 17;
APD50_iso = 4.8489; APD50_iso_std = 2.363258406; APD50_iso_sem = 0.573174355;
APD90_iso = 51.15604706; APD90_iso_std = 23.24883969; APD90_iso_sem = 5.638671864;
delta_Ca_iso = 0.463223529; delta_Ca_iso_std = 0.321491459; delta_Ca_iso_sem = 0.077973132;
Ca_tau_iso = 123.1154; Ca_tau_iso_std = 52.32411918; Ca_tau_iso_sem = 12.69046295;

figure,set(gcf,'color','w')
subplot(2,3,1),hold on
bar([APD90;APD90_iso]','FaceColor',color_mouse)
errorbar([1; 2]',[APD90;APD90_iso]',[APD90_std;APD90_iso_std]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD90 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});

subplot(2,3,2),hold on
bar([APD50;APD50_iso]','FaceColor',color_mouse)
errorbar([1; 2]',[APD50;APD50_iso]',[APD50_std;APD50_iso_std]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('ADP50 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});
title('Mouse (1-Hz) from Edwards et al. 2014')

subplot(2,3,4),hold on
bar([delta_Ca;delta_Ca_iso]','FaceColor',color_mouse)
errorbar([1; 2]',[delta_Ca;delta_Ca_iso]',[delta_Ca_std;delta_Ca_iso_std]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('CaTamp (a.u.)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});

subplot(2,3,5),hold on
bar([Ca_tau;Ca_tau_iso]','FaceColor',color_mouse)
errorbar([1; 2]',[Ca_tau;Ca_tau_iso]',[Ca_tau_std;Ca_tau_iso_std]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('CaTtau (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});
% *******************************************************************

% APD90 APD50 CaTamp CaTtau abs(RMP) APamp
% mouse control
mouse_control = [31.66589412 3.689505882 0.2054 259.9954471 80.59762941 124.2750235];
mouse_control_std = [15.4112023 1.890940989 0.133654606 106.1922159 4.56165292 13.27427913];
mouse_control_sem = [3.737765582 0.458620555 0.032416003 25.75539545 1.106363342 3.219485586];
% n = 17

% mouse iso
mouse_iso = [51.15604706 4.8489 0.463223529 123.1154 82.31813529 128.7838588];
mouse_iso_std = [23.24883969 2.363258406 0.321491459 52.32411918 4.098311014 12.73064571];
mouse_iso_sem = [5.638671864 0.573174355 0.077973132 12.69046295 0.993986423 3.087635113];
% n = 17

if n_features_mouse == 2
    mouse_control = mouse_control(1:2);
    mouse_control_std = mouse_control_std(1:2);
    mouse_control_sem = mouse_control_sem(1:2);
    mouse_iso = mouse_iso(1:2);
    mouse_iso_std = mouse_iso_std(1:2);
    mouse_iso_sem = mouse_iso_sem(1:2);
    rabbit_control = rabbit_control(1:2);
    rabbit_control_std = rabbit_control_std(1:2);
    rabbit_control_sem = rabbit_control_sem(1:2);
    rabbit_iso = rabbit_iso(1:2);
    rabbit_iso_std = rabbit_iso_std(1:2);
    rabbit_iso_sem = rabbit_iso_sem(1:2);
elseif n_features_mouse == 4
    mouse_control = mouse_control(1:4);
    mouse_control_std = mouse_control_std(1:4);
    mouse_control_sem = mouse_control_sem(1:4);
    mouse_iso = mouse_iso(1:4);
    mouse_iso_std = mouse_iso_std(1:4);
    mouse_iso_sem = mouse_iso_sem(1:4);
    rabbit_control = rabbit_control(1:4);
    rabbit_control_std = rabbit_control_std(1:4);
    rabbit_control_sem = rabbit_control_sem(1:4);
    rabbit_iso = rabbit_iso(1:4);
    rabbit_iso_std = rabbit_iso_std(1:4);
    rabbit_iso_sem = rabbit_iso_sem(1:4);
end

%% Experimentally-observed effect
rabbit_n = 6; mouse_n = 17;

ratio_mouse_std = (mouse_iso./mouse_control).*sqrt((mouse_iso_std./mouse_iso).^2+(mouse_control_std./mouse_control).^2);
ratio_rabbit_std = (rabbit_iso./rabbit_control).*sqrt((rabbit_iso_std./rabbit_iso).^2+(rabbit_control_std./rabbit_control).^2);
    ratio_mouse_sem = ratio_mouse_std./sqrt(mouse_n);
    ratio_rabbit_sem = ratio_rabbit_std./sqrt(rabbit_n);
    
% figure,set(gcf,'color','w')
% subplot(1,3,1),hold on
% b = bar([mouse_control;rabbit_control]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% if length(mouse_control) == 2
%     errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',[mouse_control;rabbit_control]',[mouse_control_std;rabbit_control_std]','k.')
% elseif length(mouse_control) == 4
%     errorbar([1-0.15 2-0.15 3-0.15 4-0.15; 1+0.15 2+0.15 3+0.15 4+0.15]',[mouse_control;rabbit_control]',[mouse_control_std;rabbit_control_std]','k.')
% end
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Exp - Control'),ylabel('APD (ms)')
% legend('Mouse','Rabbit')
% if length(mouse_control) == 2
%     set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
% elseif length(mouse_control) == 4 
% 	set(gca, 'XTick', 1:1:4, 'XTickLabel', {'APD90','APD50','CaTamp','CaTtau'})
% end
% %xlim([0.5 2.5])
% 
% subplot(1,3,2),hold on
% b = bar([mouse_iso;rabbit_iso]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% if length(mouse_control) == 2
%     errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',[mouse_iso;rabbit_iso]',[mouse_iso_std;rabbit_iso_std]','k.')
% elseif length(mouse_control) == 4
%     errorbar([1-0.15 2-0.15 3-0.15 4-0.15; 1+0.15 2+0.15 3+0.15 4+0.15]',[mouse_iso;rabbit_iso]',[mouse_iso_std;rabbit_iso_std]','k.')
% end
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Exp - ISO'),ylabel('APD (ms)')
% if length(mouse_control) == 2
%     set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
% elseif length(mouse_control) == 4 
% 	set(gca, 'XTick', 1:1:4, 'XTickLabel', {'APD90','APD50','CaTamp','CaTtau'})
% end
% %xlim([0.5 2.5])
% 
% %figure,set(gcf,'color','w')
% subplot(1,3,3),hold on
% b = bar(100*[mouse_iso./mouse_control;rabbit_iso./rabbit_control]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% if length(mouse_control) == 2
%     errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',100*[mouse_iso./mouse_control;rabbit_iso./rabbit_control]',100*[ratio_mouse_std;ratio_rabbit_std]','k.')
% elseif length(mouse_control) == 4
%     errorbar([1-0.15 2-0.15 3-0.15 4-0.15; 1+0.15 2+0.15 3+0.15 4+0.15]',100*[mouse_iso./mouse_control;rabbit_iso./rabbit_control]',100*[ratio_mouse_std;ratio_rabbit_std]','k.')
% end
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Changes after ISO'),ylabel('APD drug / APD basal (%)')
% if length(mouse_control) == 2
%     set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
% elseif length(mouse_control) == 4 
% 	set(gca, 'XTick', 1:1:4, 'XTickLabel', {'APD90','APD50','CaTamp','CaTtau'})
% end
% %xlim([0.5 2.5])

%% Simulated data (control)
mouse_sim = mean(good_outputs_X);
mouse_sim_std = std(good_outputs_X);
rabbit_sim = mean(good_outputs_Y);
rabbt_sim_std = std(good_outputs_Y);

% figure,set(gcf,'color','w'),hold on
% %subplot(1,3,1),hold on
% b = bar([mouse_sim;rabbit_sim]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% if length(mouse_control) == 2
%     errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',[mouse_sim;rabbit_sim]',[mouse_sim_std;rabbt_sim_std]','k.')
% elseif length(mouse_control) == 4
%     errorbar([1-0.15 2-0.15 3-0.15 4-0.15; 1+0.15 2+0.15 3+0.15 4+0.15]',[mouse_sim;rabbit_sim]',[mouse_sim_std;rabbt_sim_std]','k.')
% end
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Model - Control'),ylabel('APD (ms)')
% legend('Mouse','Rabbit')
% if length(mouse_control) == 2
%     set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
% elseif length(mouse_control) == 4 
% 	set(gca, 'XTick', 1:1:4, 'XTickLabel', {'APD90','APD50','CaTamp','CaTtau'})
% end
% %xlim([0.5 2.5])

%% Prediction
% Normalized ISO-effect
actual_outputs = rabbit_iso./rabbit_control.*(mean(good_outputs_Y))%;
actual_inputs = mouse_iso./mouse_control.*(mean(good_outputs_X))%;
ratio_rabbit_std_scaled = ratio_rabbit_std.*(mean(good_outputs_Y));
ratio_mouse_std_scaled = ratio_mouse_std.*(mean(good_outputs_X));

inputs = actual_inputs;
x = log(inputs);
xz = (x-mean(X))./std(X);

yz = xz*Bpls;
%yz = xz;
y = yz.*std(Y)+mean(Y);
predicted_outputs = exp(y)%;

figure,set(gcf,'color','w'),hold on
b = bar([actual_inputs;predicted_outputs;actual_outputs]','FaceColor',[1 1 1]);
set(b(1),'FaceColor',color_mouse); set(b(3),'FaceColor',color_rabbit);
if length(mouse_control) == 2
    errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',[actual_inputs;predicted_outputs;actual_outputs]',[ratio_mouse_std_scaled; 0 0 ;ratio_rabbit_std_scaled]','k.')
elseif length(mouse_control) == 4
    errorbar([1-0.22 2-0.22 3-0.22 4-0.22; 1 2 3 4; 1+0.22 2+0.22 3+0.22 4+0.22]',[actual_inputs;predicted_outputs;actual_outputs]',[ratio_mouse_std_scaled; 0 0 0 0;ratio_rabbit_std_scaled]','k.')
end
set(gca,'box','off','tickdir','out','fontsize',14)
legend('Mouse-ISO','Prediction','Rabbit-ISO')
title('"Scaled" values after drug')
ylabel('APD (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
xlim([0.5 2.5])

figure,set(gcf,'color','w'),hold on
b = bar(100*[actual_inputs./mean(good_outputs_X);predicted_outputs./mean(good_outputs_Y);actual_outputs./mean(good_outputs_Y)]','FaceColor',[1 1 1]);
set(b(1),'FaceColor',color_mouse); set(b(3),'FaceColor',color_rabbit);
if length(mouse_control) == 2
    errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',100*[actual_inputs./mean(good_outputs_X);predicted_outputs./mean(good_outputs_Y);actual_outputs./mean(good_outputs_Y)]',...
        100*[ratio_mouse_std; 0 0 ;ratio_rabbit_std]','k.')
elseif length(mouse_control) == 4
    errorbar([1-0.22 2-0.22 3-0.22 4-0.22; 1 2 3 4; 1+0.22 2+0.22 3+0.22 4+0.22]',100*[actual_inputs./mean(good_outputs_X);predicted_outputs./mean(good_outputs_Y);actual_outputs./mean(good_outputs_Y)]',...
    100*[ratio_mouse_std; 0 0 0 0;ratio_rabbit_std]','k.')
end
set(gca,'box','off','tickdir','out','fontsize',14)
legend('Mouse','Prediction','Rabbit')
title('Relative ISO-induced changes')
ylabel('APD iso / APD basal (%)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
xlim([0.5 2.5])

% figure,set(gcf,'color','w'),hold on
% b = bar(100*[predicted_outputs./mean(good_outputs_Y);actual_outputs./mean(good_outputs_Y)]','FaceColor',[1 1 1]);
% set(b(2),'FaceColor',color_rabbit);
% if length(mouse_control) == 2
%     errorbar([1+0.145 2+0.145]',100*[actual_outputs./mean(good_outputs_Y)]',100*[ratio_rabbit_std]','k.')
% elseif length(mouse_control) == 4
%     errorbar([1+0.145 2+0.145 3+0.145 4+0.145]',100*[actual_outputs./mean(good_outputs_Y)]',100*[ratio_rabbit_std]','k.')
% end
% set(gca,'box','off','tickdir','out','fontsize',14)
% legend('Prediction','Rabbit')
% title('Changes after ISO')
% ylabel('APD iso / APD basal (%)')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
% xlim([0.5 2.5])
