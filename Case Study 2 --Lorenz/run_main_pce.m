%% Create Flow-Map Approximation using Polynomial Chaos Expansions
clearvars 
clc
clear all

%Select how many transitions to use for training
Ndata = 500;
simulate_data(Ndata)
%% Define Problem
%Define whether a model already exists
%close all
trainedmodel = 0;
if trainedmodel==0;
fprintf('Beggining Training...')
uqlab
Nstates = 3;
Nparams = 3; 

load('data.mat');
Xdata = xdata;
Ydata = ydata;



Xo = Xdata(:,:) ;
Yp = Ydata(:,:);
Yo = Yp;
% Yo(:,1) = Yo(:,1) - Xo(:,1);
% Yo(:,2) = Yo(:,2) - Xo(:,2);
% Yo(:,3) = Yo(:,3) - Xo(:,3);

%Load Validation Data
load('validation_data.mat')

Xtraj = xdata(:, :);
Ytraj = ydata(:, :);


%Shuffle training data
%Fix RNG
rng(550) %550
rp = randperm(size(Xo, 1));
X = Xo(rp,:) ;
Y = Yo(rp,:) ;
% X = Xo;
% Y = Yo;
%Create Model
[pce,eltimelist] = pc_map(X,Y,Nstates,Nparams,Ndata);
save('pce_model.mat','pce')
end


sv = 4999;%249;
%Indices to compare the validation trajectories
indices = zeros(2,2);
%Enter validation indices 
indices(1,1)=1;indices(1,2)=sv;
indices(2,1)=1*sv+1;indices(2,2)=2*sv;
%indices(3,1)=2*sv+1;indices(3,2)=3*sv;
% indices(4,1)=3*sv+1;indices(4,2)=4*sv;
% indices(5,1)=4*sv+1;indices(5,2)=5*sv;
% 
% indices(6,1)=5*sv+1;indices(6,2)=6*sv;
% indices(7,1)=6*sv+1;indices(7,2)=7*sv;
% indices(8,1)=7*sv+1;indices(8,2)=8*sv;
% indices(9,1)=8*sv+1;indices(9,2)=9*sv;
% indices(10,1)=9*sv+1;indices(10,2)=10*sv;

close all

y1t_list = {};
y2t_list = {};
y3t_list = {};

y1p_list = {};
y2p_list = {};
y3p_list = {};

time_list = {};

for vv=1:length(indices)
ind1 = indices(vv,1); ind2=indices(vv,2);
npoints = ind2-ind1;
prev_ind = indices(vv,1)-1; %last index before the new traj
xinit = Xtraj(prev_ind+1,1:Nstates);

y1_true = Xtraj(prev_ind+1:prev_ind+1+npoints,1);
y2_true = Xtraj(prev_ind+1:prev_ind+1+npoints,2);
y3_true = Xtraj(prev_ind+1:prev_ind+1+npoints,3);

%y1_true = Ytraj(prev_ind+1,1);
%y2_true = Ytraj(prev_ind+1,2);


xs = Xtraj(prev_ind+1,(Nstates+1):(Nstates+Nparams));

y1 = zeros(npoints,1);y1(1) = xinit(1);
y2 = zeros(npoints,1);y2(1) = xinit(2);
y3 = zeros(npoints,1);y3(1) = xinit(3);


%% Compute pce traj
tic
tinit = 0;
time_arr = [];
time_arr = [time_arr, tinit];
for l=2:npoints+1;
% flow map mode
delta = 0.01;
tinit = tinit+delta;
time_arr = [time_arr, tinit];
% y1(l) = y1(l-1)+ (uq_evalModel(pce{1,1},[xinit,xs]));
% y2(l) = y2(l-1)+ (uq_evalModel(pce{2,1},[xinit,xs]));
% y3(l) = y3(l-1)+ (uq_evalModel(pce{3,1},[xinit,xs]));

y1(l) = (uq_evalModel(pce{1,1},[xinit,xs]));
y2(l) = (uq_evalModel(pce{2,1},[xinit,xs]));
y3(l) = (uq_evalModel(pce{3,1},[xinit,xs]));

xinit = [y1(l),y2(l),y3(l)];
end
time = 0 ;
pctime(vv) = toc;

%% Plot Trajectories
jetcolors=jet;


figure(1)
plot(time_arr,y1_true,'o','MarkerSize',5,'MarkerFaceColor','k')
hold on
%plot(time_arr,[xinit(1),y1_true],'o')
plot(time_arr,y1,'-','Color','r','LineWidth',2)
hold on
%plot(time_arr,[xinit(1),y1],'o')
box 'on'


figure(2)
plot(time_arr,y2_true,'o','MarkerSize',5,'MarkerFaceColor','k')
hold on
plot(time_arr,y2,'-','Color','g','LineWidth',2)
hold on
box 'on'



figure(3)
plot(time_arr,y3_true,'o','MarkerSize',5,'MarkerFaceColor','k')
hold on
plot(time_arr,y3,'-','Color','b','LineWidth',2)
hold on
box 'on'



figure(5)
plot(y1_true, y2_true, 'o','MarkerSize',5,'MarkerFaceColor','k')
hold on
%plot(time_arr,[xinit(1),y1_true],'o')
plot(y1,y2,'-','Color','r','LineWidth',2)
hold on
%plot(time_arr,[xinit(1),y1],'o')
box 'on'



figure(6)
plot(y1_true, y3_true, 'o','MarkerSize',5,'MarkerFaceColor','k')
hold on
%plot(time_arr,[xinit(1),y1_true],'o')
plot(y1,y3,'-','Color','g','LineWidth',2)
hold on
%plot(time_arr,[xinit(1),y1],'o')
box 'on'

figure(7)
plot(y2_true, y3_true, 'o','MarkerSize',5,'MarkerFaceColor','k')
hold on
%plot(time_arr,[xinit(1),y1_true],'o')
plot(y2,y3,'-','Color','b','LineWidth',2)
hold on
%plot(time_arr,[xinit(1),y1],'o')
box 'on'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1t_list{end+1} = y1_true;
y2t_list{end+1} = y2_true;
y3t_list{end+1} = y3_true;

y1p_list{end+1} = y1;
y2p_list{end+1} = y2;
y3p_list{end+1} = y3;

time_list{end+1} = time_arr


%% Calculate the mean relative error for each trajectory

%Loop over time-stamps
Nk = size(y1,1);
error_traj = 0.0;
err_curr = 0.0;
for k=1:Nk
xapp = [y1(k,1),y2(k,1),y3(k,1)];
xtrue = [y1_true(k,1),y2_true(k,1),y3_true(k,1)];
xdiff = xtrue - xapp;
err_curr = norm(xdiff)/norm(xtrue);
error_traj = error_traj + err_curr;
end
mte(vv,1) = error_traj/Nk;
end

%total_mean_rel_error(Nsamples) = mean(mte);



save('validations.mat', 'y1t_list','y2t_list','y3t_list',...
    'y1p_list','y2p_list','y3p_list','time_list','mte')

