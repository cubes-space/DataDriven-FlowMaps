%% Create Flow-Map Approximation using Polynomial Chaos Expansions
%function [y1ve,y2ve] = get_mte(Ndata,noise)

function mmte = get_mte(Ndata,noise)

simulate_data(Ndata,noise)
%% Define Problem
%Define whether a model already exists
%close all
trainedmodel = 0;
if trainedmodel==0;
fprintf('Beggining Training...')
uqlab
Nstates = 2;
Nparams = 1; 

load('data.mat');
Xdata = xdata;
Ydata = ydata;



Xo = Xdata(:,:) ;
Yp = Ydata(:,:);
Yo = Yp;
% Yo(:,1) = Yo(:,1) - Xo(:,1);
% Yo(:,2) = Yo(:,2) - Xo(:,2);

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


sv = 999;%249;
%Indices to compare the validation trajectories
indices = zeros(3,2);
%Enter validation indices 
indices(1,1)=1;indices(1,2)=sv;
indices(2,1)=1*sv+1;indices(2,2)=2*sv;
indices(3,1)=2*sv+1;indices(3,2)=3*sv;
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

y1p_list = {};
y2p_list = {};

time_list = {};

for vv=1:length(indices)
ind1 = indices(vv,1); ind2=indices(vv,2);
npoints = ind2-ind1;
prev_ind = indices(vv,1)-1; %last index before the new traj
xinit = Xtraj(prev_ind+1,1:Nstates);

y1_true = Xtraj(prev_ind+1:prev_ind+1+npoints,1);
y2_true = Xtraj(prev_ind+1:prev_ind+1+npoints,2);
%y1_true = Ytraj(prev_ind+1,1);
%y2_true = Ytraj(prev_ind+1,2);


xs = Xtraj(prev_ind+1,(Nstates+1):(Nstates+Nparams));

y1 = zeros(npoints,1);y1(1) = xinit(1);
y2 = zeros(npoints,1);y2(1) = xinit(2);
 
y1err = zeros(npoints,1); y1err(1) = 0;
y2err = zeros(npoints,1); y2err(1) = 0;

%% Compute pce traj
tic
tinit = 0;
time_arr = [];
time_arr = [time_arr, tinit];
for l=2:npoints+1;
% flow map mode
delta = 0.2;
tinit = tinit+delta;
time_arr = [time_arr, tinit];
% y1(l) = y1(l-1)+ (uq_evalModel(pce{1,1},[xinit,xs]));
% y2(l) = y2(l-1)+ (uq_evalModel(pce{2,1},[xinit,xs]));

y1(l) = (uq_evalModel(pce{1,1},[xinit,xs]));
y2(l) = (uq_evalModel(pce{2,1},[xinit,xs]));

%y1_rein(l) = (uq_evalModel(pce{1,1},[[y1_true(l-1),y2_true(l-1)],xs]));
%y2_rein(l) = (uq_evalModel(pce{2,1},[[y1_true(l-1),y2_true(l-1)],xs]));

%Record for the MSE
%y1err(l) = ((y1_rein(l) - y1_true(l))^2)/(y1_true(l))^2;
%y2err(l) = ((y2_rein(l) - y2_true(l))^2)/(y2_true(l))^2;


xinit = [y1(l),y2(l)];
end
time = 0 ;
pctime(vv) = toc;

%% Plot Trajectories
jetcolors=jet;

figure(1)
plot(time_arr,y1_true,'o','MarkerSize',5,'MarkerFaceColor','k')
%plot(time_arr,[xinit(1),y1_true],'o')
hold on
plot(time_arr,y1,'-','Color','r','LineWidth',2)
%plot(time_arr,[xinit(1),y1],'o')
hold on
box 'on'




figure(2)
plot(time_arr,y2_true,'o','MarkerSize',5,'MarkerFaceColor','k')
hold on
plot(time_arr,y2,'-','Color','g','LineWidth',2)
hold on
box 'on'


figure(3)
plot(y1_true, y2_true, 'o','MarkerSize',5,'MarkerFaceColor','k')
%plot(time_arr,[xinit(1),y1_true],'o')
hold on
plot(y1,y2,'-','Color','r','LineWidth',2)
%plot(time_arr,[xinit(1),y1],'o')
hold on
box 'on'



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1t_list{end+1} = y1_true;
y2t_list{end+1} = y2_true;

y1p_list{end+1} = y1;
y2p_list{end+1} = y2;

time_list{end+1} = time_arr

y1errT(vv,1) = mean(y1err)
y2errT(vv,1) = mean(y2err)

%% Calculate the mean relative error for each trajectory

%Loop over time-stamps
Nk = size(y1,1);
error_traj = 0.0;
err_curr = 0.0;
for k=1:Nk
xapp = [y1(k,1),y2(k,1)];
xtrue = [y1_true(k,1),y2_true(k,1)];
xdiff = xtrue - xapp;
err_curr = norm(xdiff)/norm(xtrue);
error_traj = error_traj + err_curr;
end
mte(vv,1) = error_traj/Nk;
end

%total_mean_rel_error(Nsamples) = mean(mte);

%Need to define the errors


%save('validations.mat', 'y1t_list','y2t_list',...
%    'y1p_list','y2p_list','time_list','mte')
mmte = mean(mte);

%y1ve = mean(y1errT);
%y2ve = mean(y2errT);