%% Create Flow-Map Approximation using Polynomial Chaos Expansions
function [pce, y1t_list,y2t_list,y3t_list,y1p_list,y2p_list,y3p_list,time_list,mte] = main_pce(Ntt);
%% Define Problem
%Define whether a model already exists
%close all
trainedmodel = 0;
if trainedmodel==0


uqlab

Nstates = 3;
Nparams = 8; 
%200 and 500 train works
Xdata = load('inputs.csv');
Ydata = load('outputs.csv');
save("Train_Data.mat",'Xdata','Ydata')

train_starts=[]
current_params = Xdata(1,Nstates+1:Nstates+Nparams);
train_starts = [1]
train_ends = []

for i=1:size(Xdata,1)
    train_p_test = Xdata(i,Nstates+1:Nstates+Nparams);

    if train_p_test == current_params;
        
    else
        train_starts=[train_starts,i];
        train_ends =[train_ends, i-1];
        current_params = Xdata(i,Nstates+1:Nstates+Nparams);
 
    end 
end 


Xo = Xdata(1:train_starts(Ntt),:) ;
Yp = Ydata(1:train_starts(Ntt),:);
Yo = Yp;
% Yo(:,1) = Yo(:,1) - Xo(:,1);
% Yo(:,2) = Yo(:,2) - Xo(:,2);
% Yo(:,3) = Yo(:,3) - Xo(:,3);

%Load Validation Data
load('10_val_ext.mat')

Xtraj = Xdata;
Ytraj = Ydata;


%Shuffle training data
%Fix RNG
rng('shuffle') %550
rp = randperm(size(Xo, 1));
X = Xo(rp,:) ;
Y = Yo(rp,:) ;


%Create Model
[pce,eltimelist] = pc_map(X,Y,Nstates,Nparams);
save('pce_model.mat','pce')
end

return

%Indices to compare the validation trajectories
indices = zeros(10,2);
%Enter validation indices 
indices(1,1)=1;  indices(1,2)=47;
indices(2,1)=48; indices(2,2)=96;
indices(3,1)=97; indices(3,2)=150;
indices(4,1)=151; indices(4,2)=202;
indices(5,1)=203; indices(5,2)=253;
indices(6,1)=254; indices(6,2)=306;

indices(7,1)=307; indices(7,2)=353;
indices(8,1)=354; indices(8,2)=402;
indices(9,1)=403; indices(9,2)=453;
indices(10,1)=454; indices(10,2)=500;


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
delta = Xtraj(prev_ind+l-1,12);
tinit = tinit+delta;
time_arr = [time_arr, tinit];
% y1(l) = y1(l-1)+ abs(uq_evalModel(pce{1,1},[xinit,xs,delta]));
% y2(l) = y2(l-1)+ abs(uq_evalModel(pce{2,1},[xinit,xs,delta]));
% y3(l) = y3(l-1)+ abs(uq_evalModel(pce{3,1},[xinit,xs,delta]));
% y4(l) = y4(l-1)+ abs(uq_evalModel(pce{4,1},[xinit,xs,delta]));

y1(l) = abs(uq_evalModel(pce{1,1},[xinit,xs,delta]));
y2(l) = abs(uq_evalModel(pce{2,1},[xinit,xs,delta]));
y3(l) = abs(uq_evalModel(pce{3,1},[xinit,xs,delta]));

xinit = [y1(l),y2(l),y3(l)];
end
time = 0 ;
pctime(vv) = toc;

%% Plot Trajectories
jetcolors=jet;

figure(1)
plot(time_arr,y1_true,'o','MarkerSize',5,'MarkerFaceColor','k')
hold on
plot(time_arr,y1,'-','Color','r','LineWidth',2)
hold on
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

%Need to define the errors


%save('validations.mat', 'y1t_list','y2t_list','y3t_list',...
 %    'y1p_list','y2p_list','y3p_list','time_list')

