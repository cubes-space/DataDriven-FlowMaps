function simulate_data(N)
%close all
%clearvars

n_runs = N+50;
sampling_time=0.01;
n_samps = 1;
sigma_min = 8; sigma_max = 12;
beta_min = 1.0  ; beta_max = 5.5;
rho_min = 10; rho_max = 30;


xdata = zeros( (n_runs-1)*n_samps,6);
ydata = zeros( (n_runs-1)*n_samps,3);

for i=1:n_runs


sigma =( sigma_max - sigma_min).*rand(1,1) + sigma_min; %
beta =( beta_max - beta_min).*rand(1,1) + beta_min; %
rho = ( rho_max - rho_min).*rand(1,1) + rho_min;


%if n_samps == 1;
tspan=[0,sampling_time];
%else
%tspan = linspace(0,samplin_time,n_samps);
%end

init1 =(10 - (-10)).*rand(1,1) -10;
init2 =(10 - (-10)).*rand(1,1) -10;
init3 =(10 - (-10)).*rand(1,1) -10; 

init_cond = [init1,init2,init3];


%[t,x] = ode45(@(t,y) model(t,y,theta), tspan, init_cond);
f = @(t,x) [-sigma*x(1) + sigma*x(2); rho*x(1) - x(2) - x(1)*x(3); -beta*x(3) + x(1)*x(2)];
[t,x] = ode45(f,tspan,init_cond);

% % figure(1)
% % %plot(tspan,x(:,1),'-o')
% % plot(tspan,[x(1,1),x(end,1)],'-o')

% % % hold on
% % % figure(2)
% % % %plot(tspan,x(:,2),'-o')
% % % plot(tspan,[x(1,2),x(end,2)],'-o')

% % hold on
% % figure(3)
% % plot(x(:,1),x(:,2),'-o')
% % hold on
%xdata ((i-1)*(n_samps-1)+1:(i*(n_samps-1)),1:2 ) = x(1,:) ;%x(1:n_samps-1,:)  ;
%xdata ((i-1)*(n_samps-1)+1:(i*(n_samps-1)),3 ) = Iapp  ;
%ydata((i-1)*(n_samps-1)+1:(i*(n_samps-1)),1:2 ) = x(end,:); %x(2:n_samps,:);
xdata(i,1:3) = x(1,:);
xdata(i,4:6) = [sigma,beta,rho];
ydata(i,1:3) = x(end,:);
end

save('data.mat', 'xdata', 'ydata') 


