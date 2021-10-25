%% Translate Wang et al. SNS data - Mouse-to-Rabbit

clear
close all
clc

disp('Mouse-to-Rabbit translation of SNS-induced effect')

%color_human = [0 0.45 0.74];
color_rabbit = [0.85 0.33 0.1];
color_mouse = [0.47 0.67 0.19];
color = [0 0 0]; % black

%% Outputs selection
% 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
% 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
% 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
% 16) CaSRmax 17) CaSRmin 18) CaSRamp
output_selection = [5 14]; % APD90 CaTtau

%% Load simulated datasets
load outputs_matrix_1500_0p1_300s_mouse_COMMON_4p8Hz
all_outputs_mouse_control_4p8Hz = all_outputs;

load outputs_matrix_1500_0p1_60s_ISO_mouse_COMMON_6p2Hz
all_outputs_mouse_iso_6p2Hz_60 = all_outputs;

load outputs_matrix_1500_0p1_300s_rabbit_COMMON_2Hz
all_outputs_rabbit_control_2Hz = all_outputs;

load outputs_matrix_1500_0p1_60s_ISO_rabbit_COMMON_3p5Hz
all_outputs_rabbit_iso_3p5Hz_60 = all_outputs;

% Populations selection
Nval = 400;

% Assemble matrices
disp('Mouse to Rabbit')

Xblock_ctrl = all_outputs_mouse_control_4p8Hz;
Yblock_ctrl = all_outputs_rabbit_control_2Hz;

Xblock_iso = all_outputs_mouse_iso_6p2Hz_60;
Yblock_iso = all_outputs_rabbit_iso_3p5Hz_60;

Xblock_ctrl = Xblock_ctrl(:,output_selection);
Xblock_iso = Xblock_iso(:,output_selection);
Yblock_ctrl = Yblock_ctrl(:,output_selection);
Yblock_iso = Yblock_iso(:,output_selection);

output_names = output_names(output_selection);
output_units = output_units(output_selection);

output_names_X = output_names;
output_units_X = output_units;

N_outputs_X = length(output_names_X);

output_names_Y = output_names;
output_units_Y = output_units;

% Check basic properties, and define X and Y
% Exclude cells with outputs = 0
good_trials = (Xblock_ctrl(:,1)>0).*(Yblock_ctrl(:,1)>0.*Xblock_iso(:,1)>0).*(Yblock_iso(:,1)>0);
good_trials_log = logical(good_trials);
% Good count: number of good trials
good_count_total = sum(good_trials_log);
% Good_outputs: array with parameters from good trials only
%good_outputs_X = Xblock(good_trials_log,:);
%good_outputs_Y = Yblock(good_trials_log,:);

good_outputs_X = Xblock_iso(good_trials_log,:)./Xblock_ctrl(good_trials_log,:);
good_outputs_Y = Yblock_iso(good_trials_log,:)./Yblock_ctrl(good_trials_log,:);

%X = log(good_outputs_X);
Y = log(good_outputs_Y);

%[N_trials N_outputs_X] = size(X);
[N_trials N_outputs_Y] = size(Y);

% Fitting vs Validation Groups
actual_inputs = good_outputs_X(end-Nval+1:end,:);
actual_outputs = good_outputs_Y(end-Nval+1:end,:);

good_count = good_count_total-Nval;
X = log(good_outputs_X(1:end-Nval,:));
Y = log(good_outputs_Y(1:end-Nval,:));

%% Calculate translation coefficients with linear regression
% PLS - nipals algorithm (2003)
[T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = ...
    PLS_nipals(X,Y,rank(X));

% % Calculate agreement of values predicted by regression (Yhat = Bpls*X) with original outputs (Y)
% % Assessment on log-transformed values
% SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
% SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
% R2each = SSYR./SSYT;
% avg_R2_fit = mean(R2each);
% 
% % Assessment on (normal) values
% R2ord_fit = zeros(1,N_outputs_Y);
% R2adj_fit = zeros(1,N_outputs_Y);
% rmse_fit = zeros(1,N_outputs_Y);
% for i = 1:N_outputs_Y
%     mdl = fitlm(exp(Y(:,i)),exp(Yhat(:,i)));
%     R2ord_fit(i) = mdl.Rsquared.Ordinary;
%     R2adj_fit(i) = mdl.Rsquared.Adjusted;
%     rmse_fit(i) = mdl.RMSE;
% end
% R2ord_fit;
% R2adj_fit;
% 
% R2_fit = R2adj_fit; % Values plotted in figures
% avg_R2calc_fit = mean(R2_fit);
% 
% % Residual Standard Deviation
% % Standard deviation for a normal distribution, centered on the predicted regression line,
% % representing the distribution of actually observed values
% oSD = zeros(1,N_outputs_Y);
% rSD = zeros(1,N_outputs_Y);
% for i = 1:N_outputs_Y
%     oSD(i) = std(exp(Y(:,i)));
%     rSD(i) = sqrt(sum((exp(Yhat(:,i)) - exp(Y(:,i)) ).^2) / (good_count-2));
% end

%% Experimental data

% APD80 CaTD80 CaTamp(n)
mouse_n = 5;
rabbit_n = 4;

% mouse control 4.8 Hz
mouse_control = [50.52 76.7 1];
mouse_control_std = [4.515196563 2.013703057 0];
mouse_control_sem = [4.515196563/sqrt(mouse_n) 2.013703057/sqrt(mouse_n) 0/sqrt(mouse_n)];

% mouse iso 6.2 Hz
mouse_iso = [49.54 69.5 1.066666667];
mouse_iso_std = [4.93791454 4.718580295 0.074475947];
mouse_iso_sem = [4.93791454/sqrt(mouse_n) 4.718580295/sqrt(mouse_n) 0.074475947/sqrt(mouse_n)];

% rabbit control 2 Hz
rabbit_control = [173.2 199.4 1];
rabbit_control_std = [14.99711083 18.97524703 0];
rabbit_control_sem = [14.99711083/sqrt(rabbit_n) 18.97524703/sqrt(rabbit_n) 0/sqrt(rabbit_n)];

% rabbit iso 3.5 Hz
rabbit_iso = [141.875 165.65 1.11];
rabbit_iso_std = [4.246861586 20.15564437 0.066332496];
rabbit_iso_sem = [4.246861586/sqrt(rabbit_n) 20.15564437/sqrt(rabbit_n) 0.066332496/sqrt(rabbit_n)];

% Heart rates
mouse_HR = [283.8 371.2];
mouse_HR_std = [40.86196275 48.77704378];
mouse_HR_sem = mouse_HR_std./sqrt(mouse_n);

rabbit_HR = [129 194.75];
rabbit_HR_std = [7.071067812 8.770214745];
rabbit_HR_sem = rabbit_HR_std./sqrt(rabbit_n);

% Figure
figure,set(gcf,'color','w')
subplot(2,3,1),hold on
bar([mouse_control(1);mouse_iso(1)]','FaceColor',color_mouse)
errorbar([1; 2]',[mouse_control(1);mouse_iso(1)]',[mouse_control_std(1);mouse_iso_std(1)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});
title('Mouse')

subplot(2,3,2),hold on
bar([mouse_control(2);mouse_iso(2)]','FaceColor',color_mouse)
errorbar([1; 2]',[mouse_control(2);mouse_iso(2)]',[mouse_control_std(2);mouse_iso_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('CaTD (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});

subplot(2,3,4),hold on
bar([rabbit_control(1);rabbit_iso(1)]','FaceColor',color_rabbit)
errorbar([1; 2]',[rabbit_control(1);rabbit_iso(1)]',[rabbit_control_std(1);rabbit_iso_std(1)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});
title('Rabbit')

subplot(2,3,5),hold on
bar([rabbit_control(2);rabbit_iso(2)]','FaceColor',color_rabbit)
errorbar([1; 2]',[rabbit_control(2);rabbit_iso(2)]',[rabbit_control_std(2);rabbit_iso_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('CaTD (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});

subplot(2,3,3),hold on
bar(1/60*[mouse_HR(1);mouse_HR(2)]','FaceColor',color_mouse)
errorbar([1; 2]',1/60*[mouse_HR(1);mouse_HR(2)]',1/60*[mouse_HR_std(1);mouse_HR_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('HR (Hz)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});

subplot(2,3,6),hold on
bar(1/60*[rabbit_HR(1);rabbit_HR(2)]','FaceColor',color_rabbit)
errorbar([1; 2]',1/60*[rabbit_HR(1);rabbit_HR(2)]',1/60*[rabbit_HR_std(1);rabbit_HR_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('HR (Hz)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','ISO'});

%% Calculate relative changes

% Experiments
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

ratio_mouse_std = (mouse_iso./mouse_control).*sqrt((mouse_iso_std./mouse_iso).^2+(mouse_control_std./mouse_control).^2);
ratio_rabbit_std = (rabbit_iso./rabbit_control).*sqrt((rabbit_iso_std./rabbit_iso).^2+(rabbit_control_std./rabbit_control).^2);
    ratio_mouse_sem = ratio_mouse_std./sqrt(mouse_n);
    ratio_rabbit_sem = ratio_rabbit_std./sqrt(rabbit_n);
    
% figure,set(gcf,'color','w')
% subplot(1,3,1),hold on
% b = bar([mouse_control;rabbit_control]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',[mouse_control;rabbit_control]',[mouse_control_std;rabbit_control_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Exp - Control')%,ylabel('APD (ms)')
% %ylim([0 300])
% legend('Mouse 4.8-Hz','Rabbit 2-Hz')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});
% 
% subplot(1,3,2),hold on
% b = bar([mouse_iso;rabbit_iso]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',[mouse_iso;rabbit_iso]',[mouse_iso_std;rabbit_iso_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Exp - ISO')%,ylabel('APD (ms)')
% %ylim([0 300])
% legend('Mouse 6.2-Hz','Rabbit 3.5-Hz')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});
% 
% %figure,set(gcf,'color','w')
% subplot(1,3,3),hold on
% b = bar(100*[mouse_iso./mouse_control;rabbit_iso./rabbit_control]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',100*[mouse_iso./mouse_control;rabbit_iso./rabbit_control]',100*[ratio_mouse_std;ratio_rabbit_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Changes after ISO'),ylabel('Ratio ISO/Control (%)')
% %ylim([0 450])
% legend('Mouse','Rabbit')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});

% Simulations
mouse_sim = mean(good_outputs_X);
mouse_sim_std = std(good_outputs_X);
rabbit_sim = mean(good_outputs_Y);
rabbt_sim_std = std(good_outputs_Y);

mouse_sim_iso = mean(Xblock_iso(good_trials_log,:));
mouse_sim_control = mean(Xblock_ctrl(good_trials_log,:));
rabbit_sim_iso = mean(Yblock_iso(good_trials_log,:));
rabbit_sim_control = mean(Yblock_ctrl(good_trials_log,:));
mouse_sim_iso_std = std(Xblock_iso(good_trials_log,:));
mouse_sim_control_std = std(Xblock_ctrl(good_trials_log,:));
rabbit_sim_iso_std = std(Yblock_iso(good_trials_log,:));
rabbit_sim_control_std = std(Yblock_ctrl(good_trials_log,:));

% figure,set(gcf,'color','w')
% subplot(1,3,1),hold on
% b = bar([mouse_sim_control;rabbit_sim_control]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',[mouse_sim_control;rabbit_sim_control]',[mouse_sim_control_std;rabbit_sim_control_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Model - Control')%,ylabel('APD (ms)')
% %ylim([0 300])
% legend('Mouse 4.8-Hz','Rabbit 2-Hz')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});
% 
% subplot(1,3,2),hold on
% b = bar([mouse_sim_iso;rabbit_sim_iso]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',[mouse_sim_iso;rabbit_sim_iso]',[mouse_sim_iso_std;rabbit_sim_iso_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Model - ISO')%,ylabel('APD (ms)')
% %ylim([0 300])
% legend('Mouse 6.2-Hz','Rabbit 3.5-Hz')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});
% 
% %figure,set(gcf,'color','w')
% subplot(1,3,3),hold on
% b = bar(100*[mouse_sim;rabbit_sim]');
% set(b(1),'FaceColor',color_mouse); set(b(2),'FaceColor',color_rabbit);
% errorbar([1-0.15 2-0.15; 1+0.15 2+0.15]',100*[mouse_sim;rabbit_sim]',100*[mouse_sim_std;rabbt_sim_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% title('Changes after ISO'),ylabel('Ratio ISO/Control (%)')
% %ylim([0 300])
% legend('Mouse','Rabbit')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});

%% Perform translation of relative SNS-induced effect
actual_outputs = rabbit_iso./rabbit_control%.*(mean(good_outputs_Y))%;
actual_inputs = mouse_iso./mouse_control%.*(mean(good_outputs_X))%;
ratio_rabbit_std_scaled = ratio_rabbit_std;%.*(mean(good_outputs_Y));
ratio_mouse_std_scaled = ratio_mouse_std;%.*(mean(good_outputs_X));

inputs = actual_inputs;
x = log(inputs);
xz = (x-mean(X))./std(X);

yz = xz*Bpls;
%yz = xz;
y = yz.*std(Y)+mean(Y);
predicted_outputs = exp(y)%;

figure,set(gcf,'color','w'),hold on
b = bar(100*[actual_inputs;predicted_outputs;actual_outputs]','FaceColor',[1 1 1]);
set(b(1),'FaceColor',color_mouse); set(b(3),'FaceColor',color_rabbit);
errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',100*[actual_inputs;predicted_outputs;actual_outputs]',100*[ratio_mouse_std_scaled; 0 0;ratio_rabbit_std_scaled]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
legend('Mouse','Prediction','Rabbit')
title('Changes after ISO')
ylabel('Output(ISO)/Output(ctrl) (%)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});
%xlim([0.5 2.5])

% figure,set(gcf,'color','w'),hold on
% b = bar(100*[predicted_outputs;actual_outputs]','FaceColor',[1 1 1]);
% set(b(2),'FaceColor',color_rabbit);
% errorbar([1+0.145 2+0.145]',100*[actual_outputs]',100*[ratio_rabbit_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% legend('Prediction','Rabbit')
% title('Changes after ISO')
% ylabel('Output(ISO)/Output(ctrl) (%)')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});
% %xlim([0.5 2.5])

%% Application to absolute experimental values
figure,set(gcf,'color','w'),hold on
b = bar([rabbit_control;predicted_outputs.*rabbit_control;rabbit_iso]','FaceColor',[1 1 1]);
set(b(1),'FaceColor',color); set(b(3),'FaceColor',color_rabbit);
errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',[rabbit_control;predicted_outputs.*rabbit_control;rabbit_iso]',[rabbit_control_std; 0 0; rabbit_iso_std]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
legend('Rabbit Control','Predicted Rabbit ISO','Rabbit ISO')
title('Application to Exp Rabbit Control data')
ylabel('Duration (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});
%xlim([0.5 2.5])

% % Application to absolute simulated values
% figure,set(gcf,'color','w'),hold on
% b = bar([rabbit_sim_control;predicted_outputs.*rabbit_sim_control;rabbit_sim_iso]','FaceColor',[1 1 1]);
% set(b(1),'FaceColor',color); set(b(3),'FaceColor',color_rabbit);
% errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',[rabbit_sim_control;predicted_outputs.*rabbit_sim_control;rabbit_sim_iso]',[rabbit_sim_control_std; 0 0; rabbit_sim_iso_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% legend('Rabbit Control','Predicted Rabbit ISO','Rabbit ISO')
% ylabel('Duration (ms)')
% title('Application to SIM Rabbit Control data')
% ylabel('Duration (ms)')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD','CaTD'});
% %xlim([0.5 2.5])
