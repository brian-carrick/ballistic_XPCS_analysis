clc; close all; clear all


% %%Figure 1
% XPCS_SANS_Fitting("LargeG2-averaging-P4-10wv-35C-a00010", 3,'P4-10 w/v%');     %These ones are good and match
% XPCS_SANS_Fitting("LargeG2-averaging-PC10P-10wv-35C-a0008", 5,'PC10P-10 w/v%');


params_list = [];

%% SI Figure 1
[q, I, params] = XPCS_SANS_Fitting("LargeG2-averaging-P4-10wv-35C-a00010", 1,'C_{10}(PC_{10})_4: 10 w/v%');     %These ones are good and match
params_list = [params_list; params];
[q, I, params] = XPCS_SANS_Fitting("LargeG2-averaging-P4-15wv-35C-a00008", 2,'C_{10}(PC_{10})_4: 15 w/v%');
params_list = [params_list; params];
[q, I, params] = XPCS_SANS_Fitting("LargeG2-averaging-P4-20wv-35C-a0006", 3,'C_{10}(PC_{10})_4: 20 w/v%');
params_list = [params_list; params];


[q, I, params] = XPCS_SANS_Fitting("LargeG2-averaging-PC10P-10wv-35C-a0008", 4,'PC_{10}P: 10 w/v%');
params_list = [params_list; params];
[q, I, params] = XPCS_SANS_Fitting("LargeG2-averaging-PC10P-12p5-35C-a00008", 5,'PC_{10}P: 12.5 w/v%');
params_list = [params_list; params];

%Figure 1
% XPCS_SANS_Fitting("LargeG2-averaging-PC10P-10wv-35C-a0008", 5,'PC_{10}P: 10 w/v%');
% XPCS_SANS_Fitting("LargeG2-averaging-P4-10wv-35C-a00010", 2,'C_{10}(PC_{10})_4: 10 w/v%');     %These ones are good and match





% XPCS_Fitting("LargeG2-averaging-P4-20wv-6C-a00006", 4, 1,'6 C')
% XPCS_Fitting("LargeG2-averaging-P4-20wv-15C-a00006", 4, 3,'15 C')
% XPCS_Fitting("LargeG2-averaging-P4-20wv-25C-a00006", 4, 5,'25 C')
% XPCS_Fitting("LargeG2-averaging-P4-20wv-35C-a0006", 4, 6,'35 C')
function [q,I,params] = XPCS_SANS_Fitting(filename, color_idx,name)
    %% Use Normally
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
    
    %%%% Downselected colors
    % color = [73 91 145;
    %     50 59 109;
    %    121 161 220;
    %     97 125 182;
    %     159 95 186;
    %     159 95 186;
    %     179 30 129;
    %     179 30 129;
    %     179 30 129;
    %     179 30 129]./256;
   
    
   
    
    file = "Data\";
    file = strcat(file,filename);
    
    data = readcell(strcat(file,"_saxs.csv"));
    data = cell2mat(data);
     
    %Pull out the q list on top & the tau's
    q = data(:,1);
    I = data(:,2);


    figure(1)
    hold on

    
    
    x = logspace(-3,-1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% Porod  %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if color_idx <= 3
        [params, errors] = fit_correlationlength_porod_fixed(q,I)   
    else
        [params, errors] = fit_correlationlength_porod(q,I)   
    end
    y = params(1)./x.^params(2) + params(3) ./ (1 + (x.*params(4)).^params(5));
    %xline(1/params(4),'-','Color',color(color_idx*2,:),'LineWidth',2,'HandleVisibility','off')
    

    %%% Plot intersections

    % Define the two terms
    term1 = @(x) params(1)./x.^params(2);
    term2 = @(x) params(3)./(1 + (x.*params(4)).^params(5));
    
    % Define function whose root gives the intersection point
    diff_fun = @(x) term1(x) - term2(x);
    
    % Provide an initial guess (adjust depending on parameter range)
    x0 = 0.01;
    
    % Use fzero to find where the two terms are equal
    x_intersect = fzero(diff_fun, x0);
    y_intersect = term1(x_intersect);
    %plot(x_intersect,6e-6,'o','Color',color(color_idx*2,:),'LineWidth',2,'MarkerSize',2,'HandleVisibility','off')
    

        %% Calculate scattering invariant
        Q_density = 1/(2*pi()*pi());
        Q_density = Q_density .* trapz(q, q.^2 .* params(3) ./ (1 + (q.*params(4)).^params(5)));
        
        Q_porod = 1/(2*pi()*pi());
        Q_porod = Q_porod .* trapz(q, q.^2 .* params(1)./q.^params(2));
        
        Ratio = Q_porod/(Q_porod+Q_density);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%% Exponential  %%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     name
% %     [params, errors] = fit_correlationlength_exp(q,I);   
% %     %params = [1e-3 600 1*1e-5 60 2]
% %     y = params(1)*exp(-(x.*params(2)).^params(3)) + params(4) ./ (1 + (x.*params(5)).^params(6));
% %     xline(1/params(2),'-','Color',color(color_idx*2,:),'LineWidth',2,'HandleVisibility','off')
% %     
% % 
% % %     %%% Plot intersections
% % % 
% % %     % Define the two terms
% % %     term1 = @(x) params(1)./x.^params(2);
% % %     term2 = @(x) params(3)./(1 + (x.*params(4)).^params(5));
% % %     
% % %     % Define function whose root gives the intersection point
% % %     diff_fun = @(x) term1(x) - term2(x);
% % %     
% % %     % Provide an initial guess (adjust depending on parameter range)
% % %     x0 = 0.01;
% % %     
% % %     % Use fzero to find where the two terms are equal
% % %     x_intersect = fzero(diff_fun, x0)
% % %     y_intersect = term1(x_intersect);
% % %     plot(x_intersect,6e-6,'o','Color',color(color_idx*2,:),'LineWidth',2,'MarkerSize',2,'HandleVisibility','off')
% % %     
% % % 
% % %         %% Calculate scattering invariant
% % %         Q_density = 1/(2*pi()*pi());
% % %         Q_density = Q_density .* trapz(q, q.^2 .* params(3) ./ (1 + (q.*params(4)).^params(5)));
% % %         
% % %         Q_porod = 1/(2*pi()*pi());
% % %         Q_porod = Q_porod .* trapz(q, q.^2 .* params(1)./q.^params(2));
% % %         
% % %         Ratio = Q_porod/(Q_porod+Q_density)
% % % 
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



       
%     %%%% broadpeak
%     name
%     [params, errors] = fit_broadpeak_porod(q,I)   
%     y = params(1)./x.^params(2) + params(3) ./ (1 + ((abs(x-params(4))*params(5)).^params(6)));
%         %LengthScales
%     xline(params(4),'-','Color',color(color_idx*2,:),'LineWidth',2)
%     xline(1/params(5),'--','Color',color(color_idx*2,:),'LineWidth',2)



%%%%%%%%% Vertically shifted plotting
%       plot(q*10,I*10^(5-color_idx),'o', 'LineWidth', 1.5, 'Color', color(color_idx*2,:),'MarkerFaceColor',color(color_idx*2,:), 'DisplayName', name);
%     if color_idx ~= 3
%         plot(x*10,y*10^(5-color_idx),'--','Color',[255 160 122]./256,'LineWidth',2,'HandleVisibility','off')
%     end
%     set(gca,'xscale','log','yscale','log','fontsize',15,'FontName','Times New Roman','linewidth',1.5,'ticklength',[.02 .02],'units',"centimeters", 'Position', [2.5 2.5 8 12.9])
%%% Not shifted plotting
    plot(q*10,I,'o', 'LineWidth', 1.5, 'Color', color(color_idx*2,:),'MarkerFaceColor',color(color_idx*2,:), 'DisplayName', name);
    plot(x*10,y,'--','Color',[255 160 122]./256,'LineWidth',2,'HandleVisibility','off')
    
    set(gca,'xscale','log','yscale','log','fontsize',15,'FontName','Times New Roman','linewidth',1.5,'ticklength',[.02 .02],'units',"centimeters", 'Position', [2.5 2.5 10 8])

    box('on');
    xlabel(strcat('\itq \rm(nm^{-1})'))
    ylabel('Intensity (cm^{-1})')
    legend('Location','northeast','fontsize',12)
    axis([0.02 0.8 3e-6 1e-3])
    xticks([0.03 0.1 0.3])

    plot([0.04 0.05], 0.00000003./([0.04 0.05].^3),'-k','LineWidth',2,'HandleVisibility','off')

    params = [params; errors];
end



function [params, errors] = fit_correlationlength(q,I)
    lb = [ 0 0];
    ub = [ inf 100 ];
    initial_guess = [1 12];
    
    f1 = fit(q, (I),'(A*1e-5 / (1 + (x*B)^2))', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region','MaxFunEvals',1e10,'MaxIter',1e10,'TolFun', 1e-20,'TolX',1e-15);
    params = coeffvalues(f1);
    params(1) = params(1)*1e-5;
    errors = confint(f1,0.95);
end

function [params, errors] = fit_correlationlength_porod(q,I)
    lb = [0 2 0 3 1];
    ub = [inf inf inf 100 2];
    initial_guess = [4 3  1 6 2];
    
    f1 = fit(q, log10(I),'log10(A*1e-11/x^B + C*1e-5 / (1 + (x*D)^E))', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region','MaxFunEvals',1e10,'MaxIter',1e10,'TolFun', 1e-30,'TolX',1e-15);
    params = coeffvalues(f1);
    params(1) = params(1)*1e-11;
    params(3) = params(3)*1e-5;
    errors = confint(f1,0.95);
    errors(:,1) = errors(:,1)*1e-11
    errors(:,3) = errors(:,3)*1e-5;
end

function [params, errors] = fit_correlationlength_porod_fixed(q,I)
    lb = [0 2 0 6 2];
    ub = [inf inf inf 100 2];
    initial_guess = [4 3  1 6 2];
    
    f1 = fit(q, log10(I),'log10(A*1e-11/x^B + C*1e-5 / (1 + (x*D)^E))', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region','MaxFunEvals',1e10,'MaxIter',1e10,'TolFun', 1e-30,'TolX',1e-20);
    params = coeffvalues(f1);
    params(1) = params(1)*1e-11;
    params(3) = params(3)*1e-5;
    errors = confint(f1,0.95);
    errors(:,1) = errors(:,1)*1e-11
    errors(:,3) = errors(:,3)*1e-5;
end

function [params, errors] = fit_correlationlength_exp(q,I)
    lb = [0 0 0.5 0 3 1];
    ub = [inf inf 2 inf 100 2];
    initial_guess = [1.2 635 1 1 17 2];
    
    f1 = fit(q, log10(I),'log10(A*1e-3*exp(-(x*B)^C) + D*1e-5 / (1 + (x*E)^F))', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region','MaxFunEvals',1e10,'MaxIter',1e10,'TolFun', 1e-30,'TolX',1e-15);
    params = coeffvalues(f1);
    params(1) = params(1)*1e-3;
    params(4) = params(4)*1e-5;
    errors = confint(f1,0.95);
end


function [params, errors] = fit_broadpeak_porod(q,I)
    lb = [0 0 0 1e-3 3 0];
    ub = [inf inf inf inf 100 4];
    initial_guess = [4 3  1 0.01 30 2];
    
    f1 = fit(q, log10(I),'log10(A*1e-11/x^B + C*1e-5 / (1 + ((abs(x-D)*E)^F)))', 'Lower', lb, 'Upper', ub, 'startpoint',initial_guess,'algorithm','trust-region','MaxFunEvals',1e10,'MaxIter',1e10,'TolFun', 1e-30,'TolX',1e-15);
    params = coeffvalues(f1);
    params(1) = params(1)*1e-11;
    params(3) = params(3)*1e-5;
    errors = confint(f1,0.95);
    errors(:,1) = errors(:,1)*1e-11
    errors(:,3) = errors(:,3)*1e-5;
end