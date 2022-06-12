function simulate_data(N,noise)
%close all
%clearvars
Cm  = 20;
gl  = 2;
Vl  =-60;
gca = 4;
Vca = 120;
gk  = 8;
Vk  = -84;
V1  = -1.2;
V2  = 18;
V3  = 12;
V4  = 17.4;
phi = 0.066;

Iapp_min = 0;
Iapp_max = 300;
n_runs = N+50;
sampling_time=0.2;
n_samps = 1;

xdata = zeros( (n_runs-1)*n_samps,3);
ydata = zeros( (n_runs-1)*n_samps,2);

for i=1:n_runs

Iapp =( Iapp_max - Iapp_min).*rand(1,1) + Iapp_min; %


theta(1) = Cm ;
theta(2) = gl;
theta(3) = Vl; 
theta(4) = gca;
theta(5) = Vca;
theta(6) = gk ;
theta(7) = Vk; 
theta(8) = V1;
theta(9) = V2;
theta(10) = V3;
theta(11) = V4;
theta(12) = phi;
theta(13) = Iapp;
%if n_samps == 1;
tspan=[0,sampling_time];
%else
%tspan = linspace(0,samplin_time,n_samps);
%end

init1 =( 75 - (-75)).*rand(1,1)  - 75; %
init2 =( 1 - 0.0).*rand(1,1)  + 0.0; %

init_cond = [init1,init2];


[t,x] = ode45(@(t,y) model(t,y,theta), tspan, init_cond);

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
xdata(i,1:2) = x(1,:);
xdata(i,3) = Iapp;
ydata(i,1:2) = x(end,:);

%sigma = noise*ydata(i,1:2);
e1 = normrnd(0,noise(1));
e2 = normrnd(0,noise(2));

ydata(i,1) = ydata(i,1) + e1;
ydata(i,2) = ydata(i,2) + e2;


end

save('data.mat', 'xdata', 'ydata') 


