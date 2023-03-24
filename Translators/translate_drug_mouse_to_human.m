%% Mouse-to-Human translation of drug effect on APD90 and APD50

close all
clear
clc

%% Build translator
disp('Mouse-to-Human translation of drug effect on APD90 and APD50')

translator_MtoH_2 % APD90, APD50

%% Experimental data
% Mouse - Portero et al. 2017
n_mouse = 5;
mouse_control = [98.9 4.3]; % control mouse  
mouse_control_sem = [8.4 0.7];
    mouse_control_std = mouse_control_sem*sqrt(n_mouse);
mouse_drug = [67.7 4]; % drug mouse
mouse_drug_sem = [4.9 0.7];
    mouse_drug_std = mouse_drug_sem*sqrt(n_mouse);

% Human - Ferrantini et al. 2018
n_human = 9;
human_control = [378.5 234];      % control human  
human_control_sem = [44.7 30];
    human_control_std = human_control_sem*sqrt(n_human);
human_drug = [354 221.8]; % drug human
human_drug_sem = [36.6 20];
    human_drug_std = human_drug_sem*sqrt(n_human);
    
%% Plot experimental data
color_human = [0 0.45 0.74];
%color_rabbit = [0.85 0.33 0.1];
color_mouse = [0.47 0.67 0.19];

figure,set(gcf,'color','w')
subplot(2,3,1),hold on
bar([human_control(1);human_drug(1)]','FaceColor',color_human)
errorbar([1; 2]',[human_control(1);human_drug(1)]',[human_control_std(1);human_drug_std(1)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD90 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','Drug'});
title('Human')

subplot(2,3,2),hold on
bar([human_control(2);human_drug(2)]','FaceColor',color_human)
errorbar([1; 2]',[human_control(2);human_drug(2)]',[human_control_std(2);human_drug_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD50 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','Drug'});

subplot(2,3,4),hold on
bar([mouse_control(1);mouse_drug(1)]','FaceColor',color_mouse)
errorbar([1; 2]',[mouse_control(1);mouse_drug(1)]',[mouse_control_std(1);mouse_drug_std(1)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD90 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','Drug'});
title('Mouse')

subplot(2,3,5),hold on
bar([mouse_control(2);mouse_drug(2)]','FaceColor',color_mouse)
errorbar([1; 2]',[mouse_control(2);mouse_drug(2)]',[mouse_control_std(2);mouse_drug_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD50 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','Drug'});

%% Prediction
% Experimental data
ratio_human_std = (human_drug./human_control).*sqrt((human_drug_std./human_drug).^2+(human_control_std./human_control).^2);
ratio_mouse_std = (mouse_drug./mouse_control).*sqrt((mouse_drug_std./mouse_drug).^2+(mouse_control_std./mouse_control).^2);
    ratio_human_sem = ratio_human_std./sqrt(n_human);
    ratio_mouse_sem = ratio_mouse_std./sqrt(n_mouse);
    
% Simulated data (control)
mouse_sim = mean(good_outputs_X);
human_sim = mean(good_outputs_Y);
mouse_sim_std = std(good_outputs_X);
human_sim_std = std(good_outputs_Y);

% Normalized drug-effect
actual_outputs = human_drug./human_control.*(mean(good_outputs_Y))%;
actual_inputs = mouse_drug./mouse_control.*(mean(good_outputs_X))%;
human_drug_std_scaled = ratio_human_std.*(mean(good_outputs_Y));
mouse_drug_std_scaled = ratio_mouse_std.*(mean(good_outputs_X));
    human_drug_sem_scaled = human_drug_std_scaled./sqrt(n_human);
    mouse_drug_sem_scaled = mouse_drug_std_scaled./sqrt(n_mouse);

inputs = actual_inputs;
x = log(inputs);
xz = (x-mean(X))./std(X);

yz = xz*Bpls;
y = yz.*std(Y)+mean(Y);
predicted_outputs = exp(y)%;

figure,set(gcf,'color','w'),hold on
b = bar([actual_inputs;predicted_outputs;actual_outputs]','FaceColor',[1 1 1]);
set(b(1),'FaceColor',color_mouse); set(b(3),'FaceColor',color_human);
errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',[actual_inputs;predicted_outputs;actual_outputs]',[mouse_drug_std_scaled; 0 0 ;human_drug_std_scaled]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
legend('Mouse','Human (predicted)','Human (actual)')
title('"Scaled" APD values after INaL block')
ylabel('APD (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});

% figure,set(gcf,'color','w'),hold on
% b = bar(100*[actual_inputs./mean(good_outputs_X);predicted_outputs./mean(good_outputs_Y);actual_outputs./mean(good_outputs_Y)]','FaceColor',[1 1 1]);
% set(b(1),'FaceColor',color_mouse); set(b(3),'FaceColor',color_human);
% errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',100*[actual_inputs./mean(good_outputs_X);predicted_outputs./mean(good_outputs_Y);actual_outputs./mean(good_outputs_Y)]',...
%     100*[ratio_mouse_std; 0 0 ;ratio_human_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% legend('Mouse','Human (predicted)','Human (actual)')
% title('APD changes with INaL block')
% ylabel('APD drug / APD basal (%)')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
