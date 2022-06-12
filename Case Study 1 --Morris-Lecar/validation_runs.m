close all
clearvars
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
n_runs = 3;
final_time = 200;
sampling = 0.2;
tspan = 0:sampling:final_time;
n_samps = length(tspan)-1;


xdata = zeros( (n_runs-1)*n_samps,3);
ydata = zeros( (n_runs-1)*n_samps,2);
Iapp = [0,60,150];

for i=1:n_runs

%Iapp =( Iapp_max - Iapp_min).*rand(1,1) + Iapp_min; %


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
theta(13) = Iapp(i);


%init1 =( 30 - (-30)).*rand(1,1)  - 30; %
%init2 =( 0.5 - 0.05).*rand(1,1)  + 0.05; %
init1 =-25;%( 10 - (-10)).*rand(1,1)  - 10; %
init2 =0.07;%( 0.25 - 0.075).*rand(1,1)  + 0.075; %


%init1 = -10;
%init2 = 0.1;

init_cond = [init1,init2];


[t,x] = ode45(@(t,y) model(t,y,theta), tspan, init_cond);

figure(1)
plot(tspan,x(:,1),'-o')
%plot(tspan,[x(1,1),x(end,1)],'-o')

hold on
figure(2)
plot(tspan,x(:,2),'-o')
%plot(tspan,[x(1,2),x(end,2)],'-o')

hold on
figure(3)
plot(x(:,1),x(:,2),'-o')
hold on
xdata ((i-1)*(n_samps-1)+1:(i*(n_samps-1)),1:2 ) =x(1:n_samps-1,:)  ;
xdata ((i-1)*(n_samps-1)+1:(i*(n_samps-1)),3 ) = Iapp(i)  ;
ydata((i-1)*(n_samps-1)+1:(i*(n_samps-1)),1:2 ) = x(2:n_samps,:);
%xdata(i,1:2) = x(1,:);
%xdata(i,3) = Iapp;
%ydata(i,1:2) = x(end,:);
end

save('validation_data.mat', 'xdata', 'ydata') 


