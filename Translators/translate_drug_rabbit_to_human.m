%% Rabbit-to-Human translation of drug effect on APD90 and APD50

close all
clear
clc

%% Build translator
disp('Rabbit-to-Human translation of drug effect on APD90 and APD50')

translator_RtoH_2 % APD90, APD50

%% Selection of block to test
% Change block_index to select the desired block
block_index = 1;
% 1) Translation of INaL block with GS-967
% 2) Translation of IKr block with E-4031
% 3) Translation of IKr block with Sotalol
% 4) Translation of IK1 block with BaCl2/PA-6
% 5) Translation of ICaL block with Verapamil

%% Experimental data

if block_index == 1
    disp('Translation of INaL block with GS-967')

    % Rabbit - Belardinelli et al. 2013
    n_rabbit = 4;
    rabbit_control = [184.3 147.2]; % control rabbit      
    rabbit_control_sem = [13.1 13.6];
        rabbit_control_std = rabbit_control_sem*sqrt(n_rabbit);
    rabbit_drug = [172.3 138.1]; % drug rabbit (1 uM)
    rabbit_drug_sem = [13.1 11.1]; % drug rabbit (1 uM)
        rabbit_drug_std = rabbit_drug_sem*sqrt(n_rabbit);

    % Human - Ferrantini et al. 2018
    n_human = 9;
    human_control = [378.5 234];      % control human  
    human_control_sem = [44.7 30];
        human_control_std = human_control_sem*sqrt(n_human);
    human_drug = [354 221.8]; % drug human
    human_drug_sem = [36.6 20];
        human_drug_std = human_drug_sem*sqrt(n_human);
elseif block_index == 2
    disp('Translation of IKr block with E-4031')

    % Human - Bussek et al. (values from Tveito et al.)
    n_human = 5;
    human_control = [356 255]; % control human
    human_control_sem = [20 15];
        human_control_std = human_control_sem*sqrt(n_human);
    human_drug = human_control.*[1.312 1.277]; % drug human relative changes
    human_drug_sem = human_control_sem.*[1.312 1.277];%[27.6 20.1];
        human_drug_std = human_drug_sem*sqrt(n_human);

    % Rabbit - Hegyi et al.
    n_rabbit = 9;
    rabbit_control = [221.5581439 170.1345514]; % control rabbit      
    rabbit_control_sem = [8.001801947 7.371768043];
        rabbit_control_std = rabbit_control_sem*sqrt(n_rabbit);
    rabbit_drug = [260.10124 203.49804]; % drug rabbit
    rabbit_drug_sem = [3.748317632 6.71191419]; % drug rabbit
        rabbit_drug_std = rabbit_drug_sem*sqrt(n_rabbit); 
elseif block_index == 3
    disp('Translation of IKr block with Sotalol')

    % Human - Orvos et al.
    n_human = 6;
    human_control = [301.8 233]; % control human
    human_control_sem = [19.7 17.4];
        human_control_std = human_control_sem*sqrt(n_human);
    human_drug = [387 281]; % drug human
    human_drug_sem = [27.6 20.1];
        human_drug_std = human_drug_sem*sqrt(n_human);

    % Rabbit - Orvos et al.
    n_rabbit = 5;
    rabbit_control = [156.9 116.6]; % control rabbit
    rabbit_control_sem = [10.7 13.1];
        rabbit_control_std = rabbit_control_sem*sqrt(n_rabbit);
    rabbit_drug = [244.3 181.3]; % drug rabbit
    rabbit_drug_sem = [15.9 21.5];
        rabbit_drug_std = rabbit_drug_sem*sqrt(n_rabbit);
elseif block_index == 4
    disp('Translation of IK1 block with BaCl2/PA-6')
    
    % Human - Jost et al. (BaCl2)
    n_human = 11;
    human_control = [280 280*52/71]; % control human
    human_control_sem = [12 12*52/71];
        human_control_std = human_control_sem*sqrt(n_human);
    human_drug = [294 294*57/80];
    human_drug_sem = [8 8*57/80];
        human_drug_std = human_drug_sem*sqrt(n_human);
    
    % Rabbit - Hegyi et al. (PA-6)
    n_rabbit = 10;
    rabbit_control = [221.5581439 170.1345514]; % control rabbit      
    rabbit_control_sem = [8.001801947 7.371768043];
        rabbit_control_std = rabbit_control_sem*sqrt(n_rabbit);
    rabbit_drug = [257.21142 179.798764]; % drug rabbit
    rabbit_drug_sem = [12.92680855 14.34862108]; % drug rabbit
        rabbit_drug_std = rabbit_drug_sem*sqrt(n_rabbit);
elseif block_index == 5
    disp('Translation of ICaL block with Verapamil')        
        
    % Human - Orvos et al.
    n_human = 3;
    human_control = [262.1 194.4]; % control human
    human_control_sem = [17.2 16.5];
        human_control_std = human_control_sem*sqrt(n_human);
    human_drug = [270.7 199.9]; % drug human
    human_drug_sem = [11.2 12.4];
        human_drug_std = human_drug_sem*sqrt(n_human);

    % Rabbit - Orvos et al.
    n_rabbit = 5;
    rabbit_control = [176.6 139.5]; % control rabbit
    rabbit_control_sem = [11.5 13.1];
        rabbit_control_std = rabbit_control_sem*sqrt(n_rabbit);
    rabbit_drug = [177.3 137.2]; % drug rabbit
    rabbit_drug_sem = [12.6 14];
        rabbit_drug_std = rabbit_drug_sem*sqrt(n_rabbit);
end

%% Plot experimental data
color_human = [0 0.45 0.74];
color_rabbit = [0.85 0.33 0.1];
%color_mouse = [0.47 0.67 0.19];

figure,set(gcf,'color','w')
subplot(2,3,1),hold on
bar([rabbit_control(1);rabbit_drug(1)]','FaceColor',color_rabbit)
errorbar([1; 2]',[rabbit_control(1);rabbit_drug(1)]',[rabbit_control_std(1);rabbit_drug_std(1)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD90 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','Drug'});
title('Rabbit')

subplot(2,3,2),hold on
bar([rabbit_control(2);rabbit_drug(2)]','FaceColor',color_rabbit)
errorbar([1; 2]',[rabbit_control(2);rabbit_drug(2)]',[rabbit_control_std(2);rabbit_drug_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD50 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','Drug'});

subplot(2,3,4),hold on
bar([human_control(1);human_drug(1)]','FaceColor',color_human)
errorbar([1; 2]',[human_control(1);human_drug(1)]',[human_control_std(1);human_drug_std(1)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD90 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','Drug'});
title('Human')

subplot(2,3,5),hold on
bar([human_control(2);human_drug(2)]','FaceColor',color_human)
errorbar([1; 2]',[human_control(2);human_drug(2)]',[human_control_std(2);human_drug_std(2)]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
ylabel('APD50 (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'Ctrl','Drug'});

%% Prediction
% Experimental data
ratio_rabbit_std = (rabbit_drug./rabbit_control).*sqrt((rabbit_drug_std./rabbit_drug).^2+(rabbit_control_std./rabbit_control).^2);
ratio_human_std = (human_drug./human_control).*sqrt((human_drug_std./human_drug).^2+(human_control_std./human_control).^2);
    ratio_rabbit_sem = ratio_rabbit_std./sqrt(n_rabbit);
    ratio_human_sem = ratio_human_std./sqrt(n_human);
    
% Simulated data (control)
rabbit_sim = mean(good_outputs_X);
human_sim = mean(good_outputs_Y);
rabbit_sim_std = std(good_outputs_X);
human_sim_std = std(good_outputs_Y);

% Normalized drug-effect
actual_outputs = human_drug./human_control.*(mean(good_outputs_Y))%;
actual_inputs = rabbit_drug./rabbit_control.*(mean(good_outputs_X))%;
human_drug_std_scaled = ratio_human_std.*(mean(good_outputs_Y));
rabbit_drug_std_scaled = ratio_rabbit_std.*(mean(good_outputs_X));
    human_drug_sem_scaled = human_drug_std_scaled./sqrt(n_human);
    rabbit_drug_sem_scaled = rabbit_drug_std_scaled./sqrt(n_rabbit);

inputs = actual_inputs;
x = log(inputs);
xz = (x-mean(X))./std(X);

yz = xz*Bpls;
y = yz.*std(Y)+mean(Y);
predicted_outputs = exp(y)%;

figure,set(gcf,'color','w'),hold on
b = bar([actual_inputs;predicted_outputs;actual_outputs]','FaceColor',[1 1 1]);
set(b(1),'FaceColor',color_rabbit); set(b(3),'FaceColor',color_human);
errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',[actual_inputs;predicted_outputs;actual_outputs]',[rabbit_drug_std_scaled; 0 0 ;human_drug_std_scaled]','k.')
set(gca,'box','off','tickdir','out','fontsize',14)
legend('Rabbit','Human (predicted)','Human (actual)')
title('"Scaled" APD values after current block')
ylabel('APD (ms)')
set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});

% figure,set(gcf,'color','w'),hold on
% b = bar(100*[actual_inputs./mean(good_outputs_X);predicted_outputs./mean(good_outputs_Y);actual_outputs./mean(good_outputs_Y)]','FaceColor',[1 1 1]);
% set(b(1),'FaceColor',color_rabbit); set(b(3),'FaceColor',color_human);
% errorbar([1-0.22 2-0.22; 1 2; 1+0.22 2+0.22]',100*[actual_inputs./mean(good_outputs_X);predicted_outputs./mean(good_outputs_Y);actual_outputs./mean(good_outputs_Y)]',...
%     100*[ratio_rabbit_std; 0 0 ;ratio_human_std]','k.')
% set(gca,'box','off','tickdir','out','fontsize',14)
% legend('Rabbit','Human (predicted)','Human (actual)')
% title('APD changes with current block')
% ylabel('APD drug / APD basal (%)')
% set(gca, 'XTick', 1:1:2, 'XTickLabel', {'APD90','APD50'});
