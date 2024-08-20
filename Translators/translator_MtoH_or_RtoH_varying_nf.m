%% Analaysis of Mouse-to-Human and Rabbit-to-Human translation
% when reducing the number of features considered from 10 to 2.

close all
clear
clc

%% Loading simulated dataset (outputs)
load outputs_current_block_human % outputs_human_block
load outputs_current_block_rabbit % outputs_rabbit_block
load outputs_current_block_mouse % outputs_mouse_block

output_names_original = {'UV', 'APpeak', '|MDP|', 'APamp', 'APD90',...
    'APD70', 'APD50', 'APD30', 'CaTmax', 'CaTmin',...
    'CaTamp', 'CaTttp', 'CaTt50', 'CaTtau', 'Na',...
    'CaSRmax', 'CaSRmin', 'CaSRamp'};

% output_units = {'mV/ms', 'mV', 'mV', 'mV', 'ms',...
%     'ms', 'ms', 'ms', 'mM', 'mM',...
%     'mM', 'ms', 'ms', 'ms', 'mM',...
%     'mM', 'mM', 'mM'};

%% Cycle
% Colors
color_human = [0 0.45 0.74];
color_rabbit = [0.85 0.33 0.1];
color_mouse = [0.47 0.67 0.19];

output_index_array = [1 2 3 4];

all_R2_output_mouse = zeros(length(output_index_array),10);
all_R2_output_rabbit = zeros(length(output_index_array),10);

all_R2_output_index_mouse = 1;
all_R2_output_index_rabbit = 2;

for ii = 1:length(output_index_array)
    output_index = output_index_array(ii);
    
    % Output selection
    % 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
    % 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
    % 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
    % 16) CaSRmax 17) CaSRmin 18) CaSRamp
    
    if output_index == 0 % ALL
        output_selection = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
    % ALL OUTPUTS
    elseif output_index == 1
        output_selection = [1 3 4 5 7 10 11 12 13 14]; % 10
    % UV	MDP	APamp	APD90	APD50	CaTmin	CaTamp	CaTttp	CaTt50	CaTtau
    elseif output_index == 2
        output_selection = [5 7 10 11 13 14]; % 6
    % APD90	APD50	CaTmin	CaTamp	CaTt50	CaTtau
    elseif output_index == 3
        output_selection = [5 7 13 14]; % 4
    % APD90	APD50	CaTt50	CaTtau
    elseif output_index == 4
        output_selection = [5 7]; % 2
    % APD90	APD50
    end
    
    X = ['# features ',num2str(length(output_selection))];
    disp(X)

    % Populations selection: Mouse-to-Human, Rabbit-to-Human
    pop_index = [all_R2_output_index_mouse all_R2_output_index_rabbit];

    n_pop_fitting = zeros(1,length(pop_index));
    avg_R2_fitting = zeros(1,length(pop_index));
    avg_R2_validation = zeros(1,length(pop_index));  

    for j = 1:length(pop_index)
        %j
        translator_index = pop_index(j);
        
        % Loading population-based simulated data
        load outputs_matrix_1500_0p1_300s_human_COMMON
        all_outputs_human_control = all_outputs;

        load outputs_matrix_1500_0p1_300s_rabbit_COMMON
        all_outputs_rabbit_control = all_outputs;

        load outputs_matrix_1500_0p1_300s_mouse_COMMON
        all_outputs_mouse_control = all_outputs;
        
        switch translator_index
            case 1
                %disp('Mouse-to-Human')
                Xblock = all_outputs_mouse_control;
                Yblock = all_outputs_human_control;
            case 2
                %disp('Rabbit-to-Human')
                Xblock = all_outputs_rabbit_control;
                Yblock = all_outputs_human_control;
        end

        % Select outputs
        Xblock = Xblock(:,output_selection);
        Yblock = Yblock(:,output_selection);
        
        output_names = output_names(output_selection);
        output_units = output_units(output_selection);
        
        % N_outputs = length(output_names);
        
        output_names_X = output_names;
        output_units_X = output_units;
        
        N_outputs_X = length(output_names_X);
        
        output_names_Y = output_names;
        output_units_Y = output_units;
        
        % Check basic properties, and define X and Y
        % Exclude cells with outputs = 0
        good_trials = (Xblock(:,1)>0).*(Yblock(:,1)>0);
        good_trials_log = logical(good_trials);
        % Good count: number of good trials
        good_count_total = sum(good_trials_log);
        % Good_outputs: array with parameters from good trials only
        good_outputs_X = Xblock(good_trials_log,:);
        good_outputs_Y = Yblock(good_trials_log,:);
        %X = log(good_outputs_X);
        Y = log(good_outputs_Y);
        
        % % No check
        % good_count = N_trials;
        % X = log(Xblock);
        % Y = log(Yblock);
        
        %[N_trials N_outputs_X] = size(X);
        [N_trials N_outputs_Y] = size(Y);
        
        % Fitting vs Validation Groups
        Nval = 400;
        actual_inputs = good_outputs_X(end-Nval+1:end,:);
        actual_outputs = good_outputs_Y(end-Nval+1:end,:);
        
        good_count = good_count_total-Nval;
        X = log(good_outputs_X(1:end-Nval,:));
        Y = log(good_outputs_Y(1:end-Nval,:));

        % Plotting options
        color = [0 0 0]; % black

        plot_matrix = 0;
        plot_fitting = 0;
        plot_validation = 1;

        % Call the PLS routine 
        N_fitting = good_count;

        % PLS - nipals algorithm (2003)
        [T,P,W,Wstar,U,B,C,Bpls,Bpls_star,Xhat,Yhat,R2x,R2y] = ...
            PLS_nipals(X,Y,rank(X));

        % N_pars, N_outputs, N_trials, good_count

        % Goodness of fit - R^2
        % Fraction of the variance in the dependent variable which is explained by the model

        % Calculate agreement of values predicted by regression (Yhat = Bpls*X) with original outputs (Y)
        % Assessment on log-transformed values
        SSYT = sum((Y-ones(good_count,1)*mean(Y)).^2);
        SSYR = sum((Yhat-ones(good_count,1)*mean(Y)).^2);
        R2each = SSYR./SSYT;
        avg_R2_fit = mean(R2each);

        % Assessment on (normal) values
        R2ord_fit = zeros(1,N_outputs_Y);
        R2adj_fit = zeros(1,N_outputs_Y);
        rmse_fit = zeros(1,N_outputs_Y);
        for i = 1:N_outputs_Y
            mdl = fitlm(exp(Y(:,i)),exp(Yhat(:,i)));
            R2ord_fit(i) = mdl.Rsquared.Ordinary;
            R2adj_fit(i) = mdl.Rsquared.Adjusted;
            rmse_fit(i) = mdl.RMSE;
        end
        R2ord_fit;
        R2adj_fit;

        R2_fit = R2adj_fit; % Values plotted in figures
        avg_R2calc_fit = mean(R2_fit);

        % Residual Standard Deviation
        % Standard deviation for a normal distribution, centered on the predicted regression line,
        % representing the distribution of actually observed values
        oSD = zeros(1,N_outputs_Y);
        rSD = zeros(1,N_outputs_Y);
        for i = 1:N_outputs_Y
            oSD(i) = std(exp(Y(:,i)));
            rSD(i) = sqrt(sum((exp(Yhat(:,i)) - exp(Y(:,i)) ).^2) / (good_count-2));
        end

        % Plot regression coefficients
        if plot_matrix == 1
            figure; set(gcf,'color','w')
            imagesc(Bpls'); colormap jet;
            set(gca,'box','off','tickdir','out','fontsize',14)
            set(gca,'YDir','normal')
            title('Regression coefficients');
            xlabel('Outputs 1');
            ylabel('Outputs 2');
            set(gca,'XTick',(1:N_outputs_X))
            set(gca,'XTickLabel',output_names_X)
            set(gca,'YTickLabel',output_names_Y)
            set(gca,'YTick',(1:N_outputs_Y))
            rotateXLabels(gca(), 90)
            colorbar
        end

        % Scatter plot
        N_figures = ceil(N_outputs_Y/10);

        if plot_fitting == 1
            dex1 = 1;
            for figdex = 1:N_figures
                figure
                set(gcf,'color','w','Position',[50,100,1500,750])
                for i = 1:10
                    if dex1 <= N_outputs_Y
                        subplot(2,5,i),hold on
                        % Plot data points
                        factor = 1;
                        if (output_index == 1) && (i == 6 || i == 7), factor = 1e6; end
                        if (output_index == 2) && (i == 3 || i == 4), factor = 1e6; end
                        plot(factor*exp(Y(:,dex1)),factor*exp(Yhat(:,dex1)),'Marker','o','LineStyle','none','Color',color)
                        xlabel(['Actual ',output_names_Y{dex1}]) 
                        ylabel(['Predicted ',output_names_Y{dex1}])
                        title(['R^2 = ',num2str(R2_fit(dex1),4)])
                        set(gca,'box','off','tickdir','out','fontsize',12)
                        % Plot identity line
                        ylim_ind = get(gca,'ylim') ;
                        xlim_ind = get(gca,'xlim') ;
                        minpoint = min([ylim_ind(1),xlim_ind(1)]);
                        maxpoint = max([ylim_ind(2),xlim_ind(2)]);
                        plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
                        xlim([minpoint, maxpoint])
                        ylim([minpoint, maxpoint])
                        dex1 = dex1+1;
                    end
                end
            end
        end

        % Validation
        N_validation = Nval;

        predicted_outputs = actual_outputs*0;
        for i = 1:Nval
            cell_index = i;
            inputs = actual_inputs(cell_index,:);
            x = log(inputs);
            xz = (x-mean(X))./std(X);

            yz = xz*Bpls;
            %    yz = xz;
            y = yz.*std(Y)+mean(Y);
            predicted_outputs(i,:) = exp(y);
        end

        % with R^2 = 0.75, the model explains approximately 75% of the variability in the predicted variable
        R2ord = zeros(1,N_outputs_Y);
        R2adj = zeros(1,N_outputs_Y);
        rmse_val = zeros(1,N_outputs_Y);
        for i = 1:N_outputs_Y
            mdl = fitlm(actual_outputs(:,i),predicted_outputs(:,i));
            R2ord(i) = mdl.Rsquared.Ordinary;
            R2adj(i) = mdl.Rsquared.Adjusted;
            rmse_val(i) = mdl.RMSE;
        end
        R2ord;
        R2adj;
        R2 = R2adj;
        avg_R2_val = mean(R2);

        % Residual Standard Deviation
        % Standard deviation for a normal distribution, centered on the predicted regression line,
        % representing the distribution of actually observed values
        oSD_val = zeros(1,N_outputs_Y);
        rSD_val = zeros(1,N_outputs_Y);
        for i = 1:N_outputs_Y
            oSD_val(i) = std(actual_outputs(:,i));
            rSD_val(i) = sqrt(sum((predicted_outputs(:,i) - actual_outputs(:,i) ).^2) / (N_validation-2));
        end

        if plot_validation == 1
            dex1 = 1;
            for figdex = 1:N_figures
                figure
                set(gcf,'color','w','Position',[50,100,1500,750])
                for i = 1:10
                    if dex1 <= N_outputs_Y
                        subplot(2,5,i),hold on
                        % Plot data points
                        factor = 1;
                        if (output_index == 1) && (i == 6 || i == 7), factor = 1e6; end
                        if (output_index == 2) && (i == 3 || i == 4), factor = 1e6; end
                        plot(factor*actual_outputs(:,dex1),factor*predicted_outputs(:,dex1),'Marker','o','LineStyle','none')%,'Color',color)
                        xlabel(['Actual ',output_names_Y{dex1}]) 
                        ylabel(['Predicted ',output_names_Y{dex1}])
                        title(['R^2 = ',num2str(R2(dex1),4)])
                        set(gca,'box','off','tickdir','out','fontsize',14)
                        % Plot identity line
                        ylim_ind = get(gca,'ylim') ;
                        xlim_ind = get(gca,'xlim') ;
                        minpoint = min([ylim_ind(1),xlim_ind(1)]);
                        maxpoint = max([ylim_ind(2),xlim_ind(2)]);
                        plot([minpoint, maxpoint],[minpoint,maxpoint],'--k')
                        xlim([minpoint, maxpoint])
                        ylim([minpoint, maxpoint])
                        dex1 = dex1+1;
                    end
                end
            end
        end

        n_pop_fitting(j) = N_fitting;
        avg_R2_fitting(j) = avg_R2_fit;
        avg_R2_validation(j) = avg_R2_val;

        if translator_index == all_R2_output_index_mouse
            all_R2_output_mouse(ii,1:length(R2)) = R2;
        end
        if translator_index == all_R2_output_index_rabbit
            all_R2_output_rabbit(ii,1:length(R2)) = R2;
        end
    end
end

%% Collecting R2 output - data
n_output = [10 6 4 2];

% Mouse
average_R2_mouse = sum(all_R2_output_mouse')./n_output;

std10 = std(all_R2_output_mouse(1,1:10));
std6 = std(all_R2_output_mouse(2,1:6));
std4 = std(all_R2_output_mouse(3,1:4));
std2 = std(all_R2_output_mouse(4,1:2));

R2_APD90_mouse = [all_R2_output_mouse(1,4) all_R2_output_mouse(2,1) all_R2_output_mouse(3,1) all_R2_output_mouse(4,1)];
R2_APD50_mouse = [all_R2_output_mouse(1,5) all_R2_output_mouse(2,2) all_R2_output_mouse(3,2) all_R2_output_mouse(4,2)];
R2_CaTt50_mouse = [all_R2_output_mouse(1,9) all_R2_output_mouse(2,5) all_R2_output_mouse(3,3) 0];
R2_CaTtau_mouse = [all_R2_output_mouse(1,10) all_R2_output_mouse(2,6) all_R2_output_mouse(3,4) 0];
R2_CaTamp_mouse = [all_R2_output_mouse(1,7) all_R2_output_mouse(2,4) 0 0];

% Rabbit
average_R2_rabbit = sum(all_R2_output_rabbit')./n_output;

std10 = std(all_R2_output_rabbit(1,1:10));
std6 = std(all_R2_output_rabbit(2,1:6));
std4 = std(all_R2_output_rabbit(3,1:4));
std2 = std(all_R2_output_rabbit(4,1:2));

R2_APD90_rabbit = [all_R2_output_rabbit(1,4) all_R2_output_rabbit(2,1) all_R2_output_rabbit(3,1) all_R2_output_rabbit(4,1)];
R2_APD50_rabbit = [all_R2_output_rabbit(1,5) all_R2_output_rabbit(2,2) all_R2_output_rabbit(3,2) all_R2_output_rabbit(4,2)];
R2_CaTt50_rabbit = [all_R2_output_rabbit(1,9) all_R2_output_rabbit(2,5) all_R2_output_rabbit(3,3) 0];
R2_CaTtau_rabbit = [all_R2_output_rabbit(1,10) all_R2_output_rabbit(2,6) all_R2_output_rabbit(3,4) 0];
R2_CaTamp_rabbit = [all_R2_output_rabbit(1,7) all_R2_output_rabbit(2,4) 0 0];

%% Show average R2 vs decreasing nf
x_flip = (1:length(n_output));
figure(101),set(gcf,'color','w'),hold on
ylabel('R2 values (-)')
xlabel('n features')
set(gca,'box','off','tickdir','out','fontsize',14)
bar(x_flip,average_R2_mouse)

plot(x_flip,R2_APD90_mouse,'--')
plot(x_flip,R2_APD50_mouse,'--')
legend('Average','APD90','APD50')

errorbar(x_flip,average_R2_mouse,[std10 std6 std4 std2],'k.');

for jj = 1:length(n_output)
    y = all_R2_output_mouse(jj,1:n_output(jj));%x = n_output(jj)*ones(1, + rand(1,200);  
    mean(y);
    
    x = x_flip(jj)*ones(1,n_output(jj)) + 0.2*rand(1,n_output(jj)) - 0.1; 
    plot(x,y,'ko','LineWidth',1.5)%,'MarkerSize',10)
end
ylim([0 1])
xticks([1 2 3 4])
xticklabels({'n = 10','n = 6','n = 4','n = 2'})
title('Mouse-to-Human')

figure(102),set(gcf,'color','w'),hold on
ylabel('R2 values (-)')
xlabel('n features')
set(gca,'box','off','tickdir','out','fontsize',14)
bar(x_flip,average_R2_rabbit)

plot(x_flip,R2_APD90_rabbit,'--')
plot(x_flip,R2_APD50_rabbit,'--')
legend('Average','APD90','APD50')

errorbar(x_flip,average_R2_rabbit,[std10 std6 std4 std2],'k.');

for jj = 1:length(n_output)
    y = all_R2_output_rabbit(jj,1:n_output(jj));%x = n_output(jj)*ones(1, + rand(1,200);  
    mean(y);
    
    x = x_flip(jj)*ones(1,n_output(jj)) + 0.2*rand(1,n_output(jj)) - 0.1; 
    plot(x,y,'ko','LineWidth',1.5)%,'MarkerSize',10)
end
ylim([0 1])
%xlim([0.5 4.5])
xticks([1 2 3 4])
xticklabels({'n = 10','n = 6','n = 4','n = 2'})
title('Rabbit-to-Human')
