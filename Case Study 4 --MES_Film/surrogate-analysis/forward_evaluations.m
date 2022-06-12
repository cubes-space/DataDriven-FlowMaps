clearvars
clear all
clc

maxNumCompThreads('automatic');
Nrel = 20000;
uqlab
%Load the PCE model
load('pce_model.mat')
load('train_data_os.mat')
load('time_list.mat')

Xdata = X;
Ydata = Y;
%save("Train_Data.mat",'Xdata','Ydata')

samples_avail=0;

if samples_avail == 0

Xm=[]
npoints = length(Tlist);
rng('shuffle')

%Assign the integration lags
nparams=6;
timelags=zeros(length(Tlist)-1,1);
for kk=1:length(timelags)
timelags(kk) = Tlist(kk+1)-Tlist(kk);
end


%Initialize arrays
pctime = zeros(Nrel,1);
xs_rel = zeros(Nrel,nparams);
y1 = zeros(npoints,1,Nrel);
y2 = zeros(npoints,1,Nrel);
y3 = zeros(npoints,1,Nrel);
 
for nn=1:Nrel
nn
xs = zeros(nparams,1);
nstates=3;
x0=  [7.5e-5,9e-3,7e-3];
for j=1:nparams
xs(j) =(max(Xdata(:,j+nstates))-min(Xdata(:,j+nstates))).*rand(1,1) + min(Xdata(:,j+nstates)); %eg %0.309   
end

y1(1,1,nn) = x0(1);
y2(1,1,nn) = x0(2);
y3(1,1,nn) = x0(3);

%compute pce traj
tinit = 0;
time_arr_rel = [];
time_arr_rel = [time_arr_rel, tinit];

tic
for ll=2:npoints;
% flow map mode
delta = timelags(ll-1);
tinit = tinit+delta;
time_arr_rel = [time_arr_rel, tinit];
y1(ll,1,nn) = abs(uq_evalModel(pce{1,1},[x0,xs',delta]));
y2(ll,1,nn) = abs(uq_evalModel(pce{2,1},[x0,xs',delta]));
y3(ll,1,nn) = abs(uq_evalModel(pce{3,1},[x0,xs',delta]));

x0 = [y1(ll,1,nn),y2(ll,1,nn),y3(ll,1,nn)];
end

elapsed_time_toc = toc
pctime(nn) = elapsed_time_toc;

xs_rel(nn,:) = xs;


end


save('saved_reals.mat', 'Tlist', 'y1', 'y2', 'y3', 'xs_rel', 'pctime');

end
