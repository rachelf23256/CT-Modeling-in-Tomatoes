%% CT MODEL Project %%
% Revised: October 8, 2022
% Author: 

%% Deterministic model

close all
clear

numsims = 10000; % number of simulations 

t = [8, 20, 30, 33, 42, 51, 63, 70]; % time in days
y=[0, 1, 3, 6, 9, 9, 15, 21]; % Disease counts from bio data


p0 = [0.0005, 0.0001, 5]; % Initial parameter vector: P(1)-Beta beets, p(2)-Beta tom, p(3)=i.c. beets
Parameter_vector = fminsearch(@(p) norm(y-fit_beets(p, t)), p0);      % Estimate Parameters

detmodout=fit_beets(abs(Parameter_vector),t); % Output of ODE model at optimied parameters


%% Stochastic Model
Parameter_vector(3)=round(Parameter_vector(3)); % Rounding initial number of infected

% yr: output of stochastic model representing number of infected tomatoes
% tr: time steps that stochastic model outputs are for
[yr, tr]=stochmod(abs(Parameter_vector), numsims,t); % output of stochastic model 
stochmodout=mean(yr); 

%Generates figure with all three model outputs
figure(1)
plot(t, y,  'k*','linewidth', 1); 
hold on 
plot(t, detmodout, 'b+','linewidth', 1); 
hold on 
plot(t, stochmodout, 'ro','linewidth', 1)
title('Number of Infected Tomatoes Over Time'); 
xlabel('Days After Planting')
ylabel('Number of Infected Tomatoes')
legend('Experimental Data', 'Deterministic Model', 'Mean of the Stochastic Model', 'location', 'northwest'); 

 
% Generate figure for stochastic boxplot
figure(2)
scatter(t,y,'m','linewidth', 1)
hold on
boxplot(yr, t, 'symbol', '', 'Positions',t)
title('Number of Infected Tomatoes Over Time')
xlabel('Days After Planting')
ylabel('Number of Infected Tomatoes')
legend('Experimental Data', 'location', 'northwest')

%%  Evaluating Variance


t2 = 8:1:70; % Additional time values to have stochastic model output at
[yr2, tr2] = stochmod(abs(Parameter_vector), numsims,t2);
baseline_var = sqrt([var(yr2)]); % sstoring value for later use

% Plot of the standard deviation
subplot(1,2,2)
plot(t2, sqrt([var(yr2)]), 'b*');
xlabel('Days After Planting'); 
ylabel('Standard Deviation'); 
title('Standard Deviation vs. Time');

%% Sample paths
r= randi(numsims,10,1); % Random selection of paths

% Plotting paths 
subplot(1,2,1)
for h=1:10
    plot(t2,yr2(r(h), :) )
    hold on
    xlim([0,70])
end
hold on; 
title('Sample Paths of Stochastic Model'); 
xlabel('Days After Planting');
ylabel('Number of Infected Tomatoes');
%% Effects of perturbation of disease transmission parameters on mean of stochastic model
%
% Perturbing beet to beet transmission rate
original_parameters=Parameter_vector; % Holding the original parameter values  

Perturbed_b2b=zeros(8, 20); 
Perturbed_beta_prime=zeros(1,20);

% Generation of different beet to beet transmission rates
% Ranging from -50% to 50% by 5% increments
for r=1:10
   Perturbed_beta_prime(r)=-r*0.05*original_parameters(1)+original_parameters(1);
end

f=Perturbed_beta_prime(1:10); 
f=flip(f); 
Perturbed_beta_prime(1:10)=f; 

temp2=zeros(1,10);
for r=1:10
    temp2(r)=r*0.05*original_parameters(1)+original_parameters(1);
end
Perturbed_beta_prime(11:20)=temp2; 

% stochastic model output for different beet to beet transmission rates
 for r=1:20
     Parameter_vector(1)=Perturbed_beta_prime(r); 
    [yr, tr]=stochmod(abs(Parameter_vector), numsims, t);
    Perturbed_b2b(:, r)=mean(yr); 
 end


 
% Plot of mean of stochastic model with different b2b transmision rate 
figure(6)
subplot(1,2,1)
plot(t, Perturbed_b2b(:, 14), 'o--', 'linewidth', 1)
hold on
hold on 
plot(t, stochmodout, "*k", 'linewidth', 1)
plot(t, Perturbed_b2b(:,7), 'd--', 'linewidth', 1)
hold on 
ylim([0,25]);
legend('+20%', 'Baseline','-20%', 'location', 'northwest', 'FontSize',14) 
title({'Disease Counts with Perturbed', 'Beet-to-Beet Transmission Rate'})
xlabel('Days after Planting')
ylabel('Number of Infected Tomatoes')
ax = gca;
ax.FontSize = 14;
ax.ColorOrder = [0.4940 0.1840 0.5560; 0 0 0;  0.4940 0.1840 0.5560 ];




%% Perturbing beet-to-tomato tranmission rate

Parameter_vector = original_parameters; 

Perturbed_b2t = zeros(8, 20);
Perturbed_beta=zeros(1,20);

% Generation of different beet to beet transmission rates
% Ranging from -50% to 50% by 5% increments
for r=1:10
    Perturbed_beta(r)=-r*0.05*original_parameters(2)+original_parameters(2);
end

f=Perturbed_beta(1:10); 
f=flip(f); 
Perturbed_beta(1:10)=f; 

temp2=zeros(1,10);
for r=1:10
    temp2(r)=r*0.05*original_parameters(2)+original_parameters(2);
end
Perturbed_beta(11:20)=temp2; 

% stochastic model output for different b2t transmission rates
 for r=1:20
     Parameter_vector(2)=Perturbed_beta(r); 
    [yr, tr] = stochmod(abs(Parameter_vector), numsims, t);
    Perturbed_b2t(:, r)=mean(yr); 
 end

% Plot of Perturbed Beet-to-tomato transmission rate 
subplot(1,2,2)
plot(t, Perturbed_b2t(:, 14),'+--','linewidth', 1) 
hold on
plot(t, Perturbed_b2t(:, 12), '+--', 'linewidth', 1)
hold on
plot(t, Perturbed_b2t(:, 11), '+--', 'linewidth', 1)
hold on
plot(t, stochmodout, 'k*-', 'linewidth', 1)
plot(t, Perturbed_b2t(:, 10), 'o--', 'linewidth', 1)
hold on 
plot(t, Perturbed_b2t(:,9),'o--', 'linewidth', 1)
hold on 
plot(t, Perturbed_b2t(:, 7),'o--', 'linewidth', 1)
hold on 
legend('+20%', '+10%','+5%','Baseline', '-5%', '-10%','-20%', 'location', 'northwest', 'FontSize',14) 
title({'Disease Counts with perturbed', 'beet-to-tomato transmission rate'})
xlabel('Days after Planting')
ylabel('Number of Infected Tomatoes')
ax = gca;
ax.ColorOrder = [0.4940 0.1840 0.5560; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0 0 0; 0.8500 0.3250 0.0980;0 0.4470 0.7410; 0.4940 0.1840 0.5560 ];
ax.FontSize = 14;
%% Stochastic model STD and variance with perturbed beet to beet transmission rate

t2 = 8:1:70;
perturbed_b2b_var = zeros(length(t2), 20);
parameter_vector = original_parameters; 
for r=1:20
     parameter_vector(1)=Perturbed_beta_prime(r); 
    [yr, tr, variance_vals]= stochmod(abs(parameter_vector), numsims, t2);
    perturbed_b2b_var(:, r)=var(yr); 
 end

% Plot of Perturbed Beet-to-tomato transmission rate 
figure(9)
subplot(1,2,1)
plot(t2,  sqrt(perturbed_b2b_var(:, 14)), '+--', 'linewidth', 1)
hold on
plot(t2, baseline_var, 'k*', 'linewidth', 1)
hold on 
plot(t2,  sqrt(perturbed_b2b_var(:, 7)),'o--', 'linewidth', 1)
legend('-20%','Baseline', '20%', 'location', 'northwest', 'FontSize',14) 
title({'Standard deviation of stochastic model', 'with perturbed beet-to-beet', 'transmission rate'})
xlabel('Days after Planting')
ylabel('Number of Infected Tomatoes')
ax = gca;
ax.FontSize = 14;
ax.ColorOrder = [0.4940 0.1840 0.5560; 0 0 0;  0.4940 0.1840 0.5560 ];
%% Stochastic Model STD and variance with perturbed beet to tomato transmission rate

t2 = 8:1:70;
perturbed_b2t_var = zeros(length(t2), 20);
parameter_vector = original_parameters; 
for r=1:20
     parameter_vector(2)=Perturbed_beta(r); 
    [yr, tr]=stochmod(abs(parameter_vector), numsims, t2);
    perturbed_b2t_var(:, r)=var(yr); 
end

% Plot of STD of stochastic model with perturbed Beet-to-tomato transmission rate 
subplot(1,2,2)
plot(t2, sqrt(perturbed_b2t_var(:, 14)), '+--', 'linewidth', 1)
hold on
plot(t2,  sqrt(perturbed_b2t_var(:, 12)), '+--', 'linewidth', 1)
hold on
plot(t2,  sqrt(perturbed_b2t_var(:, 11)), '+--', 'linewidth', 1)
plot(t2, baseline_var, 'k*', 'linewidth', 1)
plot(t2,  sqrt(perturbed_b2t_var(:, 10)), 'o--', 'linewidth', 1)
hold on 
plot(t2,  sqrt(perturbed_b2t_var(:,9)),'o--', 'linewidth', 1)
hold on 
plot(t2,  sqrt(perturbed_b2t_var(:, 7)),'o--', 'linewidth', 1)
hold on 
legend('-20%', '-10%','-5%','Baseline', '5%', '10%','20%', 'location', 'northwest', 'FontSize',14) 
title({'Standard deviation of stochastic model', 'with perturbed beet-to-tomato', 'transmission rate'})
xlabel('Days after Planting')
ylabel('Number of Infected Tomatoes')
ax = gca;
ax.FontSize = 14;
ax.ColorOrder = [0.4940 0.1840 0.5560; 0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0 0 0; 0.8500 0.3250 0.0980;0 0.4470 0.7410; 0.4940 0.1840 0.5560 ];
