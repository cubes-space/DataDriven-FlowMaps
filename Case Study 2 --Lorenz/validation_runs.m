close all
clearvars

n_runs = 2;
final_time = 50;
sampling = 0.01;
tspan = 0:sampling:final_time;
n_samps = length(tspan)-1;

sigma_min = 8; sigma_max = 12;
beta_min = 1.0  ; beta_max = 5.5;
rho_min = 10; rho_max = 30;


xdata = zeros( (n_runs-1)*n_samps,6);
ydata = zeros( (n_runs-1)*n_samps,3);

for i=1:n_runs

sigma =( sigma_max - sigma_min).*rand(1,1) + sigma_min; %
beta =( beta_max - beta_min).*rand(1,1) + beta_min; %
rho = ( rho_max - rho_min).*rand(1,1) + rho_min;

if i==1
sigma = 10 ;
beta = 8/3;
rho = 28;
end



if i==2
sigma = 10 ;
beta = 8/3;
rho = 15;
end
%theta(1) =sigma ;
%theta(2) = rho;
%theta(3) = beta; 


init1 =(2 - (-2)).*rand(1,1) -2;
init2 =(2 - (-2)).*rand(1,1) -2;
init3 =(2 - (-2)).*rand(1,1) -2; 

%init1 = -10;
%init2 = 0.1;

init_cond = [init1,init2,init3];

f = @(t,x) [-sigma*x(1) + sigma*x(2); rho*x(1) - x(2) - x(1)*x(3); -beta*x(3) + x(1)*x(2)];
[t,x] = ode45(f,tspan,init_cond);     % Runge-Kutta 4th/5th order ODE solver
figure(1)
plot(x(:,1),x(:,2))
hold on
figure(2)
plot(x(:,1),x(:,3))
hold on
figure(3)
plot(x(:,2),x(:,3))
hold on
figure(4)
hold on
plot(x(:,1),'o')
figure(5)
hold on
plot(x(:,2),'o')
figure(6)
hold on
plot(x(:,3),'o')
figure(7)
hold on
plot3(x(:,1),x(:,2),x(:,3),'-')



% % % [t,x] = ode45(@(t,y) model(t,y,theta), tspan, init_cond);
% % % 
% % % figure(1)
% % % plot(tspan,x(:,1),'-o')
% % % %plot(tspan,[x(1,1),x(end,1)],'-o')
% % % hold on
% % % figure(2)
% % % plot(tspan,x(:,2),'-o')
% % % %plot(tspan,[x(1,2),x(end,2)],'-o')
% % % 
% % % 
% % % hold on
% % % figure(3)
% % % plot(tspan,x(:,3),'-o')
% % % %plot(tspan,[x(1,2),x(end,2)],'-o')
% % % 
% % % hold on
% % % figure(3)
% % % plot(x(:,1),x(:,2),'-o')
% % % hold on
cp = [sigma,beta,rho];
plength = n_samps-1;
p_matr = repmat(cp,plength,1);

xdata ((i-1)*(n_samps-1)+1:(i*(n_samps-1)),1:3 ) =x(1:n_samps-1,:)  ;
xdata ((i-1)*(n_samps-1)+1:(i*(n_samps-1)),4:6) = p_matr ;
ydata((i-1)*(n_samps-1)+1:(i*(n_samps-1)),1:3 ) = x(2:n_samps,:);
cp=[]
plength=[]
p_matr = []

end

save('validation_data.mat', 'xdata', 'ydata') 


