clc; close all; clear all

f
function [KWW_list, first_cumulant_list, KWW_tau_scaling, KWW_tau_scaling_uncertainty, cage_size,q_list] = XPCS_Fitting(filename, q_max_idx, color_idx,name)
    %% Use Normally
    %Use for collapse plot
    color = [21 21 47;
        27 30 73;
        50 59 109;
        73 91 145;
        97 125 182;
        121 161 220;
        131 139 216;
        144 119 205;
        159 95 186;
        171 68 161;
        179 30 129;
        136 27 93]./256;
    

    
    N_rep = 100;     % 100 Bootstrapped replicates
    
    
    
    
    
    
    
    
    
    first_cumulant_list = [];
    KWW_list = [];
    KWW_list = [];
    reduced_color_list =[];
    
    file = "Data\";
    file = strcat(file,filename);
    
    data = readcell(strcat(file,"_g2.csv"));
    data = cell2mat(data);
    uncertainty_data = readcell(strcat(file,"_g2err.csv"));
    uncertainty_data = cell2mat(uncertainty_data);
    
    
    %Pull out the q list on top & the tau's
    q_list = data(1,2:size(data,2))';
    tau_0 = data(2:size(data,1),1);
    data = data(2:size(data,1),2:size(data,2));
    uncertainty_data = uncertainty_data(2:size(uncertainty_data,1),2:size(uncertainty_data,2));
    
    
    %% Plotting
    q_list = q_list(1:q_max_idx);
    num_per_subplot = 1;                        % how many q's per subplot
    num_subplots = ceil(q_max_idx/num_per_subplot);
    
    % choose grid close to square
    ncols = ceil(sqrt(num_subplots));
    nrows = ceil(num_subplots / ncols);
    
    figure(1+color_idx)
    tiledlayout(nrows,ncols);                      % auto-wrapped grid
    
    figure(20+color_idx)
    ax_all = axes;
    hold(ax_all,'on');
    for q_idx =  1:q_max_idx
    % Which subplot should this go in?
    
    
        % Extract data
        q = q_list(q_idx);
        tau = tau_0;
        S_qt = data(:,q_idx)-1;
        err_S_qt = uncertainty_data(:,q_idx);

        %15 for good data
        %17 for temp sweep
        tau = tau(17:end);
        S_qt = S_qt(17:end);
        err_S_qt = err_S_qt(17:end);


        % Label
        q_name = strcat('\itq = \rm',string(round(q*10,2)),"nm^{-1}");
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Single Exponential Fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cumulant_parameters, error_tau] = fit_singleexponential(tau, S_qt);
        A_cumulant = cumulant_parameters(1);
        first_cumulant_tau = cumulant_parameters(2);
        B_cumulant = cumulant_parameters(3);
        first_cumulant_list = [first_cumulant_list; first_cumulant_tau, A_cumulant, B_cumulant];
    
        % Fit curve for plotting
        tau_plot = logspace(-4,1);
        first_cumulant_plot = A_cumulant.*exp(-(tau_plot./first_cumulant_tau).^(1))+B_cumulant;
    
        % Pick marker shape by position within group of 3
        markers = {'o','s','v'};
        marker = markers{mod(q_idx-1,num_per_subplot)+1};
        c = color(mod(q_idx-1,num_per_subplot)+1,:);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% KWW Exponential Fit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [KWW_parameters, KWW_uncertainties] = fit_KWW(tau, S_qt);
        
        A_KWW = KWW_parameters(1);
        KWW_tau = KWW_parameters(2);
        Beta = KWW_parameters(3);
        B_KWW = KWW_parameters(4);
        KWW_tau_avg = (KWW_tau/Beta) * gamma(1/Beta);
        if Beta == 1.5
            KWW_tau_avg_unc = (KWW_uncertainties(2)/Beta) * gamma(1/Beta);
            KWW_uncertainties(3) = 0;
        else
            %uncertainty on averaged function
            df_dtau = (1 ./ Beta) .* gamma(1 ./ Beta);
            g_term = gamma(1 ./ Beta);
            psi_term = psi(1 ./ Beta); % digamma function
            df_dbeta = -(KWW_tau ./ (Beta.^2)) .* g_term .* (1 + (1 ./ Beta) .* psi_term);
            KWW_tau_avg_unc = sqrt((df_dtau .* KWW_uncertainties(2)).^2 + (df_dbeta .* KWW_uncertainties(3)).^2);
        end


        KWW_list = [KWW_list; A_KWW, KWW_uncertainties(1), KWW_tau, KWW_uncertainties(2), Beta, KWW_uncertainties(3), KWW_tau_avg, KWW_tau_avg_unc];
    
        % Fit curve for plotting
        tau_plot = logspace(-4,1);
        y_plot = A_KWW.*exp(-(tau_plot./KWW_tau).^(Beta))+B_KWW;
    
    
    
    
    
    
        % Pick marker shape by position within group of 3
        markers = {'o','s','v'};
        marker = markers{mod(q_idx-1,num_per_subplot)+1};
        c = color(mod(q_idx-1,num_per_subplot)+1,:);
    
        % Plot
        figure(1+color_idx)
        subplot_idx = ceil(q_idx / num_per_subplot);
        ax = nexttile(subplot_idx);
        hold(ax,'on');
        errorbar(ax, tau, S_qt, err_S_qt, marker, ...
            'LineWidth', 2, 'Color', c, 'MarkerFaceColor', c, 'DisplayName', q_name);
        plot(ax, tau_plot, y_plot, '-', ...
            'LineWidth', 2, 'Color', color(8,:), 'DisplayName','\Beta = 3/2');
        plot(ax, tau_plot, first_cumulant_plot, '--', ...
            'LineWidth', 2, 'Color', color(5,:), 'DisplayName','\Beta = 1');
        % Formatting
        set(ax, 'XScale', 'log', 'YScale', 'lin', 'FontSize', 15, 'FontName', 'Times New Roman', 'LineWidth', 1.5, 'TickLength', [.02 .02]);
        box(ax,'on');
        %ylim([0.001 1])
        %legend(ax,'Location','southeast');
    

        %%%% All traces plot
        if q_idx == 3
            figure(20)
            title(strcat(string(round(q*10,2)),"nm^{-1}"))
            hold on
%                 errorbar(tau, (S_qt-B_KWW)./A_KWW, err_S_qt./A_KWW, marker, ...
%                     'LineWidth', 1.5, 'Color', color(color_idx*2,:), 'MarkerFaceColor', color(color_idx*2,:), 'DisplayName', strcat(string(round(q*10,2)),"nm^{-1}",filename));
%                 plot(tau_plot, (y_plot-B_KWW)./A_KWW, '-', 'LineWidth', 1.5, 'Color', color(color_idx*2,:), 'HandleVisibility','off');

                errorbar(tau, (S_qt-B_KWW), err_S_qt, marker, ...
                    'LineWidth', 1.5, 'Color', color(color_idx*2,:), 'MarkerFaceColor', color(color_idx*2,:), 'DisplayName', name);
                plot(tau_plot, (y_plot-B_KWW), '-', 'LineWidth', 1.5, 'Color', color(color_idx*2,:), 'HandleVisibility','off');

            
            % --- after the for-loop ends ---
            set(gca,'xscale','log','yscale','lin','fontsize',15,'FontName','Times New Roman','linewidth',1.5,'ticklength',[.02 .02],'units',"centimeters", 'Position', [2.5 2.5 10 8])
            box('on');
            legend('Location','bestoutside');
            xlabel('\tau');
            ylabel('g_2(q,t)-1');
            legend('Location','northeast','fontsize',12)
            axis([0.01 1 0 0.4])
        end


          
    % % %     % BootStrapping for Uncertainty
    % % %     S_qt_BootStrapped = bootstrap_NSE(S_qt,err_S_qt,N_rep);
    % % %     
    % % %     first_cumulant_bootstrapped = zeros(N_rep,2);
    % % %     first_cumulant_bootstrapped(1,:) = [A, first_cumulant_tau];
    % % %     for repeat = 2:N_rep
    % % %         [first_cumulant_parameters, relative_error_tau] = fit_singleexponential(tau(first_cumulant_start_idx:first_cumulant_end_idx), S_qt_BootStrapped(first_cumulant_start_idx:first_cumulant_end_idx,repeat));
    % % %         first_cumulant_bootstrapped(repeat,:) = [first_cumulant_parameters(1) first_cumulant_parameters(2)];
    % % %     end
    % % %     first_cumulant_list = [first_cumulant_list; mean(first_cumulant_bootstrapped(:,2)), std(first_cumulant_bootstrapped(:,2))]
    % % % 
    % % % 
    % % %     %%% Parameter Histograms
    % % %     figure(3)
    % % %     ax = nexttile(q_idx);    
    % % %     
    % % %     title(ax, strcat("\itq\rm = ", num2str(round(q,3)), " ", char(197), "^{-1}"));
    % % %     
    % % %     set(ax, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 14, ...
    % % %         'FontName', 'Times New Roman', 'LineWidth', 1.5, 'TickLength', [0.04 0.04]);
    % % %     box(ax, 'on');
    % % % 
    % % % 
    % % % 
    % % %     histogram(first_cumulant_bootstrapped(:,2),logspace(log10(min(first_cumulant_bootstrapped(:,2))),log10(max(first_cumulant_bootstrapped(:,2))),100),'FaceColor',color(q_idx,:))
    
    
    
    
%         % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % % % %%%%%%%%%%%%%%%%%%%%%%%%   KWW Fit   %%%%%%%%%%%%%%%%%%%%%%%%%%
%         % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%         
%             
%             [KWW_parameters, relative_error_tau] = fit_KWW(tau, S_qt);
%             KWW_A = KWW_parameters(1);
%             KWW_tau = KWW_parameters(2);
%             KWW_beta = KWW_parameters(3);
%            
%         
%             % BootStrapping for Uncertainty
%             S_qt_BootStrapped = bootstrap_NSE(S_qt,err_S_qt,N_rep);
%             
%             KWW_bootstrapped_tau_beta = zeros(N_rep,4);
%             KWW_bootstrapped_tau_beta(1,:) = [KWW_A, KWW_tau, KWW_beta,(KWW_tau/KWW_beta) * gamma(1/KWW_beta)];
%             for repeat = 2:N_rep
%                 [KWW_parameters, relative_error_tau] = fit_KWW(tau, S_qt_BootStrapped(:,repeat));
%                 if q_idx == 6 && KWW_parameters(2) > 0.7 
%                     KWW_bootstrapped_tau_beta(repeat,:) = [KWW_parameters(1) KWW_parameters(2) KWW_parameters(3) ((KWW_parameters(2)./KWW_parameters(3)) .* gamma(1./KWW_parameters(3)))];
%                 end
%             end
%     
%             for repeat = N_rep:-1:2
%                 if KWW_bootstrapped_tau_beta(repeat,1) == 0
%                     KWW_bootstrapped_tau_beta(repeat,:) = [];
%                 end
%             end
%             %KWW_list = [KWW_list; mean(KWW_bootstrapped_tau_beta(:,3)),std(KWW_bootstrapped_tau_beta(:,3))]
%             KWW_list = [KWW_list; mean(KWW_bootstrapped_tau_beta(:,1)),std(KWW_bootstrapped_tau_beta(:,1)),mean(KWW_bootstrapped_tau_beta(:,2)),std(KWW_bootstrapped_tau_beta(:,2)),mean(KWW_bootstrapped_tau_beta(:,3)),std(KWW_bootstrapped_tau_beta(:,3)),mean(KWW_bootstrapped_tau_beta(:,4)),std(KWW_bootstrapped_tau_beta(:,4))]
%     
%     
%     
%                 %%% Parameter Histograms
%         figure(30+color_idx)
%         ax = nexttile(q_idx);    
%         
%         title(ax, strcat("\tau_{KWW}, \itq\rm = ", num2str(round(q,3)), " ", char(197), "^{-1}"));
%         
%         set(ax, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 14, ...
%             'FontName', 'Times New Roman', 'LineWidth', 1.5, 'TickLength', [0.04 0.04]);
%         box(ax, 'on');
%     
%     
%     
%         histogram(KWW_bootstrapped_tau_beta(:,3),logspace(log10(min(KWW_bootstrapped_tau_beta(:,3))),log10(max(KWW_bootstrapped_tau_beta(:,3))),100),'FaceColor',color(4,:))
%     
%         figure(40+color_idx)
%         ax = nexttile(q_idx);    
%         
%         title(ax, strcat("\beta_{KWW}, \itq\rm = ", num2str(round(q,3)), " ", char(197), "^{-1}"));
%         
%         set(ax, 'XScale', 'log', 'YScale', 'linear', 'FontSize', 14, ...
%             'FontName', 'Times New Roman', 'LineWidth', 1.5, 'TickLength', [0.04 0.04]);
%         box(ax, 'on');
%     
%     
%      histogram(KWW_bootstrapped_tau_beta(:,2),logspace(log10(min(KWW_bootstrapped_tau_beta(:,2))),log10(max(KWW_bootstrapped_tau_beta(:,2))),100),'FaceColor',color(4,:))
%     
    
    end
    
    
    
        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%   Relaxation Plots   %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    figure(104)
    
    hold on
    %KWW 
    [coeff2, error_power] = fit_linear(log(q_list), log(KWW_list(:,7)));
    KWW_tau_scaling = coeff2(1)
    KWW_tau_scaling_uncertainty = error_power
    x_fit = [q_list(1), q_list(end)];
    y_fit = exp(coeff2(2)+KWW_tau_scaling.*log(x_fit));
    
    %plot(x_fit, y_fit,'-','LineWidth',2,'Color',color(color_idx*2,:),'HandleVisibility','off')
    errorbar(q_list*10,KWW_list(:,7),KWW_list(:,8),'s','MarkerSize',7, 'MarkerFaceColor',color(color_idx*2,:),'linewidth',2,'Color',color(color_idx*2,:),'DisplayName',name)
    
    %% Exponential
    % [coeff2, error_power] = fit_linear(log(q_list), log(first_cumulant_list(:,1)));
    % scaling = coeff2(1);
    % exp_tau_scaling = coeff2(1)
    % exp_tau_scaling_uncertainty = error_power
    % x_fit = [q_list(1), q_list(end)];
    % y_fit = exp(coeff2(2)+exp_tau_scaling.*log(x_fit));
    
    % plot(x_fit, y_fit,'-','LineWidth',2,'Color',color(3,:),'HandleVisibility','off')
    % plot(q_list,first_cumulant_list(:,1),'s','MarkerSize',7, 'MarkerFaceColor',color(3,:),'linewidth',2,'Color',color(3,:),'DisplayName','Exponential')
    
    
    xlabel(strcat('\itq \rm(nm^{-1})'))
    ylabel(strcat('\it\tau \rm(s)'))
    set(gca,'xscale','log','yscale','log','fontsize',14,'FontName','Times New Roman','linewidth',1.5,'ticklength',[.02 .02],'units',"centimeters", 'Position', [2 2 10 8])
    box on
    legend('Location','Southwest')
    axis([0.026 0.26 0.045 1.1])
    xticks([0.05 0.1 0.2])
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%   Amplitude Plots   %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    figure(105)
    
    hold on
    %KWW 
    [coeffs, error] = fit_Stretched_DW(q_list(2:end).^2, KWW_list(2:end,1));
    
    x_fit = logspace(log10(q_list(1)), log10(q_list(end)));
    y_fit = coeffs(1)*exp(-(x_fit.^2./coeffs(2)).^coeffs(3));
    sqrt(3./coeffs(2))
    cage_size = [sqrt(3./coeffs(2)), (error(2)/2) * sqrt(3./coeffs(2).^3)]
    Beta = coeffs(3)
    plot(x_fit*10, y_fit./coeffs(1),'-','LineWidth',2,'Color',color(color_idx*2,:),'HandleVisibility','off')
    errorbar(q_list*10,(KWW_list(:,1))./coeffs(1),KWW_list(:,2)./coeffs(1),'s','MarkerSize',7, 'MarkerFaceColor',color(color_idx*2,:),'linewidth',2,'Color',color(color_idx*2,:),'DisplayName',name)
    


%     plot(x_fit, y_fit,'-','LineWidth',2,'Color',color(color_idx*2,:),'HandleVisibility','off')
%     errorbar(q_list,(KWW_list(:,1)),KWW_list(:,2),'s','MarkerSize',7, 'MarkerFaceColor',color(color_idx*2,:),'linewidth',2,'Color',color(color_idx*2,:),'DisplayName',name)
%     


    % %% Exponential
    % [coeff2, error_power] = fit_linear(log(q_list), log(first_cumulant_list(:,2)-first_cumulant_list(:,3)))
    % scaling = coeff2(1);
    % x_fit = [q_list(1), q_list(end)];
    % y_fit = exp(coeff2(2)+scaling.*log(x_fit));
    % 
    % plot(x_fit, y_fit,'-','LineWidth',2,'Color',color(3,:),'HandleVisibility','off')
    % plot(q_list,(first_cumulant_list(:,2)-first_cumulant_list(:,3)),'s','MarkerSize',7, 'MarkerFaceColor',color(3,:),'linewidth',2,'Color',color(3,:),'DisplayName','\Beta = 1')
    
    
    xlabel(strcat('\itq \rm(nm^{-1})'))
    ylabel(strcat('\itA/A\rm(\itq \rm= 0.05 (nm^{-1}))'))
    ylabel(strcat('\itA/A_o'))
    set(gca,'xscale','log','yscale','log','fontsize',14,'FontName','Times New Roman','linewidth',1.5,'ticklength',[.02 .02],'units',"centimeters", 'Position', [2 2 10 8])
    box on
    legend('Location','Southwest','fontsize',12)
    %axis([0.04 0.25 0.1 1.1])
    axis([0.026 0.26 0.12 1.6])
    xticks([0.05 0.1 0.2])
    yticks([0.2 0.5 1])



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%   Beta Plots   %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    figure(106)
    
    hold on
    errorbar(q_list*10,(KWW_list(:,5)),KWW_list(:,6),'s','MarkerSize',7, 'MarkerFaceColor',color(color_idx*2,:),'linewidth',2,'Color',color(color_idx*2,:),'DisplayName',name)
    xlabel(strcat('\itq \rm(nm^{-1})'))
    ylabel('KWW \Beta')
    set(gca,'xscale','log','yscale','log','fontsize',14,'FontName','Times New Roman','linewidth',1.5,'ticklength',[.02 .02],'units',"centimeters", 'Position', [2 2 10 8])
    box on
    legend('Location','Southwest')


end



function [coeff1, error_slope] = fit_linear(x,y)

lb = [-inf -inf];
ub = [inf inf];
initial_guess = [-1; 1];

f1 = fit(x, y,'a*x+b', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
coeff1 = coeffvalues(f1);
slope = coeff1(1);

conf_interval_slope = confint(f1,0.95);
error_slope = [abs(conf_interval_slope(1,1)-slope),abs(conf_interval_slope(2,1)-slope)]; % 95% confidence region

end 

%%% KWW Stretched Exponential Fit
function [coeff1, relative_error_tau] = fit_KWW(tau, S_qt)


lb = [0 1e-2 1.5 0];
ub = [2 inf 1.5 0.0];
initial_guess = [0.15;0.1; 1.5; S_qt(end)];

f1 = fit(tau, S_qt,'a*exp(-(x/b)^c+d)', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
coeff1 = coeffvalues(f1);
a = coeff1(1);
tau = coeff1(2);
beta = coeff1(3);

conf_interval_tau = confint(f1,0.95);
relative_error_tau = [max(abs(conf_interval_tau(1,1)-a),abs(conf_interval_tau(2,1)-a)),     max(abs(conf_interval_tau(1,2)-tau),abs(conf_interval_tau(2,2)-tau)),   max(abs(conf_interval_tau(1,3)-beta),abs(conf_interval_tau(2,3)-beta))]; % 95% confidence region
end 

%%% KWW Stretched Exponential Fit
function [coeff1, relative_error_tau] = fit_Stretched_DW(tau, S_qt)


lb = [0 0 1];
ub = [inf inf 1];
initial_guess = [0.15;0.1; 1];

f1 = fit(tau, S_qt,'a*exp(-(x/b)^c)', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
coeff1 = coeffvalues(f1);
a = coeff1(1);
tau = coeff1(2);
beta = coeff1(3);

conf_interval_tau = confint(f1,0.95);
relative_error_tau = [max(abs(conf_interval_tau(1,1)-a),abs(conf_interval_tau(2,1)-a)),     max(abs(conf_interval_tau(1,2)-tau),abs(conf_interval_tau(2,2)-tau)),   max(abs(conf_interval_tau(1,3)-beta),abs(conf_interval_tau(2,3)-beta))]; % 95% confidence region
end 

%%%% First Cumulant Fit
function [coeff1, error_tau] = fit_singleexponential(tau, S_qt)


lb = [0 0 0];
ub = [2 inf 2.5];
initial_guess = [0.15;0.1;S_qt(end) ];

f1 = fit(tau, (S_qt),'a*exp(-(x/b).^(1)) + c', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
%f1 = fit(tau, log(S_qt),'a - (2/3)*(x/b)', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
coeff1 = coeffvalues(f1);
a = coeff1(1);
tau = coeff1(2);

conf_interval_tau = confint(f1,0.95);
error_tau = [min(abs(conf_interval_tau(1,2)-tau),abs(conf_interval_tau(2,2)-tau)), max(abs(conf_interval_tau(1,2)-tau),abs(conf_interval_tau(2,2)-tau))]; % 95% confidence region
end 
%%%% First Cumulant Fit
function [coeff1, error_tau] = fit_singleexponential_nobkg(tau, S_qt)


lb = [0 0];
ub = [2 inf];
initial_guess = [0.15;0.01.^2];

f1 = fit(tau, (S_qt),'a*exp(-(x/b))', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
%f1 = fit(tau, log(S_qt),'a - (2/3)*(x/b)', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
coeff1 = coeffvalues(f1);
a = coeff1(1);
tau = coeff1(2);

conf_interval_tau = confint(f1,0.95);
error_tau = [min(abs(conf_interval_tau(1,2)-tau),abs(conf_interval_tau(2,2)-tau)), max(abs(conf_interval_tau(1,2)-tau),abs(conf_interval_tau(2,2)-tau))]; % 95% confidence region
end 

% % % 
% % % function [coeff1, relative_error_tau] = fit_exp_bkg(tau, S_qt)
% % % 
% % % 
% % % lb = [0 0 0 0];
% % % ub = [2 inf 1 1];
% % % initial_guess = [0.5;10; 1; 0.5];
% % % 
% % % f1 = fit(tau, S_qt,'a*exp(-(x/b)^c)+d', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
% % % coeff1 = coeffvalues(f1);
% % % a = coeff1(1);
% % % tau = coeff1(2);
% % % 
% % % conf_interval_tau = confint(f1,0.95);
% % % relative_error_tau = max(abs(conf_interval_tau(1,2)-tau)/tau,abs(conf_interval_tau(2,2)-tau)/tau); % 95% confidence region
% % % end 
% % % 
% % % function [coeff1, relative_error_tau] = fit_double_exp(tau, S_qt)
% % % 
% % % 
% % % lb = [0 0 0 0];
% % % ub = [2 inf 2 inf];
% % % initial_guess = [0.5; 1; 0.5; 20];
% % % 
% % % f1 = fit(tau, S_qt,'a*exp(-(x/b))+c*exp(-(x/d))', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
% % % coeff1 = coeffvalues(f1);
% % % a = coeff1(1);
% % % tau = coeff1(2);
% % % 
% % % conf_interval_tau = confint(f1,0.95);
% % % relative_error_tau = max(abs(conf_interval_tau(1,2)-tau)/tau,abs(conf_interval_tau(2,2)-tau)/tau); % 95% confidence region
% % % end 
% % % 
% % % 
% % % function [coeff1, error_power] = fit_powerlaw(x, y)
% % % 
% % % 
% % % lb = [-inf -10 -inf];
% % % ub = [inf 10 inf];
% % % initial_guess = [1;-1; 0];
% % % 
% % % f1 = fit(x, y,'a*x^b + c', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region');
% % % coeff1 = coeffvalues(f1);
% % % a = coeff1(1);
% % % b = coeff1(2);
% % % c = coeff1(3);
% % % 
% % % conf_interval = confint(f1,0.95)
% % % error_power = max(abs(conf_interval(1,2)-b),abs(conf_interval(2,2)-b)); % 95% confidence region
% % % end 
% % % 
% % % 
function [I_rep] = bootstrap_NSE(I,dI,N_rep)
% Inputs are (1) I [inverse cm]
%             (2) dI [inverse cm] >> this is "sigma"
%             (3) N_rep is the number of boostrap replicates to do
% Output: (1) I_rep: a matrix of size N_q x N_rep, where N_q is the number
%             of q-values
I_rep = zeros(size(I,1),N_rep);

%% Create replicates
% Populate the M = 1 slice using the original data
I_rep(:,1) = I;
% Create replicas that are N(Icut,SDcut^2). Note, MATLAB reads the Gaussian
% distribution as N(Icut,SDcut), so it uses the SD, not the variance.
% Obviously, do not create replicas if M = 1.
if N_rep > 1
    for i = 2:N_rep
        I_rep(:,i) = normrnd(I,dI,size(I));
    end
end
end


