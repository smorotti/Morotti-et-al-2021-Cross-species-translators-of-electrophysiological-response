%% Recursive feature elimination analysis of Mouse-to-Human
% or Rabbit-to-Human translation (control, 1-Hz)

close all
clear
clc

%% Selection of input species

% Change input_species_index to select the input species
% 1 for Mouse, otherwise for Rabbit
input_species_index = 1;

%% Features selection
% Input features selection (Mouse or Rabbit)
% 1) UV 2) APpeak 3) -MDP 4) APamp 5) APD90
% 6) APD70 7) APD50 8) APD30 9) CaTmax 10) CaTmin
% 11) CaTamp 12) CaTttp 13) CaTt50 14) CaTtau 15) Namin
% 16) CaSRmax 17) CaSRmin 18) CaSRamp

output_index = 0;
%input_features_selection = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]; % ALL

input_features_selection = [1 3 4 5 7 10 11 12 13 14]; % 10
% UV	MDP	APamp	APD90	APD50	CaTmin	CaTamp	CaTttp	CaTt50	CaTtau

% Output features selection (Human)
output_features_selection = [5 7 11 14]; % 4
% APD90	APD50 CaTamp CaTtau

%% Loading data
% Mouse
load outputs_matrix_1500_0p1_300s_mouse_COMMON
all_outputs_mouse_control = all_outputs;

% Rabbit
load outputs_matrix_1500_0p1_300s_rabbit_COMMON
all_outputs_rabbit_control = all_outputs;

% Human
load outputs_matrix_1500_0p1_300s_human_COMMON
all_outputs_human_control = all_outputs;

if input_species_index == 1
    disp('Mouse-to-Human')
    Xblock_input = all_outputs_mouse_control;
else
    disp('Rabbit-to-Human')
    Xblock_input = all_outputs_rabbit_control;
end

Yblock_input = all_outputs_human_control;

%% Select outputs
Xblock = Xblock_input(:,input_features_selection);
input_features_names = output_names(input_features_selection);
input_features_units = output_units(input_features_selection);

% N_outputs = length(output_names);
output_names_X = input_features_names;
output_units_X = input_features_units;

N_outputs_X = length(input_features_names);

Yblock = Yblock_input(:,output_features_selection);
output_features_names = output_names(output_features_selection);
output_features_units = output_units(output_features_selection);

output_names_Y = output_features_names;
output_units_Y = output_features_units;

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

%% Regression - Whole input dataset
disp('Whole input dataset')

% Plotting options
color = [0 0 0]; % black

plot_matrix = 1;
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

n_pop_fitting = N_fitting;
avg_R2_fitting = avg_R2_fit;
avg_R2_validation = avg_R2_val;

% Final outputs
all_R2_start = R2
avg_R2_start = avg_R2_validation
std_R2_start = std(R2)

input_features_selection_start = input_features_selection;
    
%% Recursive evaluation (with elimination of least informative feature)
disp('Recursive analysis')

plot_index = 0;

% Final outputs
input_features_discard_order = zeros(1,length(input_features_selection_start));
all_R2_final = zeros(length(input_features_selection_start),N_outputs_Y);
avg_R2_final = zeros(1,length(input_features_selection_start));
std_R2_final = zeros(1,length(input_features_selection_start));

for jj = 1:length(input_features_selection_start)
    input_features_discard = zeros(1,length(input_features_selection));
    input_features_selection_roi = zeros(length(input_features_selection),length(input_features_selection)-1);
    all_R2_out = zeros(length(input_features_selection),N_outputs_Y);
    avg_R2_out = zeros(1,N_outputs_Y);
    
    for ii = 1:length(input_features_selection)
        input_features_discard(ii) = input_features_selection(ii);
        if ii == 1
            input_features_selection_roi(ii,:) = input_features_selection(2:end);
        elseif ii == length(input_features_selection)
            input_features_selection_roi(ii,:) = input_features_selection(1:end-1);
        else
            input_features_selection_roi(ii,:) = [input_features_selection(1:ii-1) input_features_selection(ii+1:end)];
        end

        % Update X
        Xblock = Xblock_input(:,input_features_selection_roi(ii,:));
        input_features_names = output_names(input_features_selection_roi(ii,:));
        input_features_units = output_units(input_features_selection_roi(ii,:));
        output_names_X = input_features_names;
        output_units_X = input_features_units;
        N_outputs_X = length(input_features_names);

        good_outputs_X = Xblock(good_trials_log,:);
        actual_inputs = good_outputs_X(end-Nval+1:end,:);
        X = log(good_outputs_X(1:end-Nval,:));
        
        % Regression
        plot_matrix = 0;
        plot_fitting = 0;
        plot_validation = 0;

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

        n_pop_fitting = N_fitting;
        avg_R2_fitting = avg_R2_fit;
        avg_R2_validation = avg_R2_val;

        all_R2_out(ii,:) = R2;
        avg_R2_out(ii) = avg_R2_validation;
    end
    
    all_R2_out;
    avg_R2_out;
    [val, index] = max(avg_R2_out);
    
    all_R2_final(jj,:) = all_R2_out(index,:);
    avg_R2_final(jj) = val;
    std_R2_final(jj) = std(all_R2_out(index,:));
    
    input_features_discard_order(jj) = input_features_discard(index);
    
    input_features_selection = input_features_selection_roi(index,:);
end

all_R2_final
avg_R2_final
std_R2_final

% Input features discard order
discard_list = output_names(input_features_discard_order);
print_discard_list = discard_list'

%% Plot results
% X-axis
iteration = (0:1:length(input_features_selection_start));
% Combined
combined_all_R2 = [all_R2_start; all_R2_final];
combined_avg_R2 = [avg_R2_start avg_R2_final];
combined_std_R2 = [std_R2_start std_R2_final];

figure,set(gcf,'color','w'),hold on
errorbar(iteration,combined_avg_R2,combined_std_R2)
ylabel('R2 values (-)')
xlabel('# iteration')
plot(iteration,combined_all_R2(:,1),'--')
plot(iteration,combined_all_R2(:,2),'--')
plot(iteration,combined_all_R2(:,3),'--')
plot(iteration,combined_all_R2(:,4),'--')
set(gca,'box','off','tickdir','out','fontsize',14)
ylim([0 1])
xlim([0 length(input_features_selection_start)-1])
legend('average',output_features_names{1},output_features_names{2},output_features_names{3},output_features_names{4})

if input_species_index == 1
    title('Mouse-to-Human');
else
    title('Rabbit-to-Human');
end
