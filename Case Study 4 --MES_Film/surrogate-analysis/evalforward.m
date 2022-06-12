function [y1,y2,y3]=evalforward(xs_rel,Nrel,Tlist,npoints)

uqlab
%Load the PCE model
load('pce_model.mat')
load('train_data_os.mat')
Xdata = X;
Ydata = Y;
xs=xs_rel';

%Assign the integration lags
timelags=zeros(length(Tlist)-1,1);
for kk=1:length(timelags)
timelags(kk) = Tlist(kk+1)-Tlist(kk);
end

%Initialize arrays
pctime = zeros(Nrel,1);
%xis_rel = zeros(Nrel,11);
y1 = zeros(npoints,1);
y2 = zeros(npoints,1);
y3 = zeros(npoints,1);

nstates=3;
x0=  [7.5e-5,9e-3,7e-3];

y1(1) = x0(1);
y2(1) = x0(2);
y3(1) = x0(3);

%compute pce traj
tinit = 0;
time_arr_rel = [];
time_arr_rel = [time_arr_rel, tinit];

tic
for ll=2:npoints
% flow map mode
delta = timelags(ll-1);
tinit = tinit+delta;
time_arr_rel = [time_arr_rel, tinit];
y1(ll) = (uq_evalModel(pce{1,1},[x0,xs',delta]));
y2(ll) = (uq_evalModel(pce{2,1},[x0,xs',delta]));
y3(ll) = (uq_evalModel(pce{3,1},[x0,xs',delta]));

x0 = [y1(ll),y2(ll),y3(ll)];
end
pctime = toc;

xs_rel = xs;


%save('saved_reals.mat', 'Tlist', 'y1', 'y2', 'y3', 'y4','xs_rel');

end


