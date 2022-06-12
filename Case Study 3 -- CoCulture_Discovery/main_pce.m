%% Create Flow-Map Approximation using Polynomial Chaos Expansions
uqlab
%% Define Problem
%Define whether a model already exists
trainedmodel = 0;   

if trainedmodel==0
clearvars
clear all
clc
close all

uqlab

Nstates = 5;
Ncontrol = 3; 
%200 and 500 train works
Xdata = load('Xtrain.csv');
Ydata = load('Ytrain.csv');

Xo = Xdata;
Yp = Ydata;
Yo = Yp;

%Load Validation Data

Xtraj = load('Xtest.csv');
Ytraj = load('Ytest.csv');


%Shuffle training data
%Fix RNG
rng('shuffle') %550
rp = randperm(size(Xo, 1));
X = Xo(rp,:) ;
Y = Yo(rp,:) ;
%Create Model
[pce,eltimelist] = pc_map(X,Y,Nstates,Ncontrol);
save('pce_model.mat','pce')
end

%Indices to compare the validation trajectories
indices = zeros(1,2);
%Enter validation indices 
indices(1,1)=1;  indices(1,2)=40;
%indices(2,1)=41; indices(2,2)=80;


% indices(3,1)=101; indices(3,2)=150;
% indices(4,1)=151; indices(4,2)=200;
% indices(5,1)=201; indices(5,2)=250;

% indices(3,1)=201; indices(3,2)=300;




close all

y1t_list = {};
y2t_list = {};
y3t_list = {};
y4t_list = {};
y5t_list = {};

y1p_list = {};
y2p_list = {};
y3p_list = {};
y4p_list = {};
y5p_list = {};

y1v_list = {};
y2v_list = {};
y3v_list = {};
y4p_list = {};
y5p_list = {};


time_list = {};

for vv=1:size(indices,1)
    ind1 = indices(vv,1); ind2=indices(vv,2);
npoints = ind2-ind1;
prev_ind = indices(vv,1)-1; %last index before the new traj
xinit = Xtraj(prev_ind+1,1:5);

y1_true = Xtraj(prev_ind+1:prev_ind+1+npoints,1);
y2_true = Xtraj(prev_ind+1:prev_ind+1+npoints,2);
y3_true = Xtraj(prev_ind+1:prev_ind+1+npoints,3);
y4_true = Xtraj(prev_ind+1:prev_ind+1+npoints,4);
y5_true = Xtraj(prev_ind+1:prev_ind+1+npoints,5);


y1= zeros(npoints,1);y1(1) = xinit(1);
y2 = zeros(npoints,1);y2(1) = xinit(2);
y3 = zeros(npoints,1);y3(1) = xinit(3);
y4 = zeros(npoints,1);y4(1) = xinit(4);
y5 = zeros(npoints,1);y5(1) = xinit(5);

y1v= zeros(npoints,1);
y2v = zeros(npoints,1);
y3v = zeros(npoints,1);
y4v = zeros(npoints,1);
y5v = zeros(npoints,1);


%% Compute pce traj
tic
tinit = 0;
time_arr = [];
inputs_seq_c1in=[]
inputs_seq_c2in=[]
inputs_seq_q = []

time_arr = [time_arr, tinit];
for l=2:npoints+1;
% flow map mode
delta = Xtraj(prev_ind+l-1,9);
tinit = tinit+delta;
time_arr = [time_arr, tinit];
xs = Xtraj(prev_ind+l,6:8);
inputs_seq_c1in = [inputs_seq_c1in,xs(1)];
inputs_seq_c2in = [inputs_seq_c2in,xs(2)];
inputs_seq_q = [inputs_seq_q,xs(3)];


% y1(l) = y1(l-1)+ abs(uq_evalModel(pce{1,1},[xinit,xs,delta]));
% y2(l) = y2(l-1)+ abs(uq_evalModel(pce{2,1},[xinit,xs,delta]));
% y3(l) = y3(l-1)+ abs(uq_evalModel(pce{3,1},[xinit,xs,delta]));
% y4(l) = y4(l-1)+ abs(uq_evalModel(pce{4,1},[xinit,xs,delta]));
% y5(l) = y5(l-1)+ abs(uq_evalModel(pce{5,1},[xinit,xs,delta]));

[y1(l),y1v(l)] = uq_evalModel(pce{1,1},[xinit,xs,delta]);
[y2(l),y2v(l)] = uq_evalModel(pce{2,1},[xinit,xs,delta]);
[y3(l),y3v(l)] = uq_evalModel(pce{3,1},[xinit,xs,delta]);
[y4(l),y4v(l)] = uq_evalModel(pce{4,1},[xinit,xs,delta]);
[y5(l),y5v(l)] = uq_evalModel(pce{5,1},[xinit,xs,delta]);
% y4(l) = abs(uq_evalModel(pce{4,1},[xinit,xs,delta]));
% y5(l) = abs(uq_evalModel(pce{5,1},[xinit,xs,delta]));

xinit = [y1(l),y2(l),y3(l),y4(l),y5(l)];


end
time = 0 ;
pctime(vv) = toc;

%% Plot Trajectories
jetcolors=jet;
hms = 3 %how many sigma 
figure(1)
y1plus = y1+hms*sqrt(y1v);
y1minus = y1-hms*sqrt(y1v);

plot(time_arr,y1_true,'s','Linewidth',2,'Color','k')
hold on
plot(time_arr,y1plus,'--','Color','r','LineWidth',2)
hold on
plot(time_arr,y1minus,'--','Color','r','LineWidth',2)
%hold on
%patch([time_arr fliplr(time_arr)], [ y1minus fliplr(y1plus)], 'r')
hold on
plot(time_arr,y1,'-','Color','r','LineWidth',2)
hold on
box 'on'

figure(2)
y2plus = y2+hms*sqrt(y2v);
y2minus = y2-hms*sqrt(y2v);

plot(time_arr,y2_true,'s','Linewidth',2,'Color','k')
hold on
plot(time_arr,y2plus,'--','Color','b','LineWidth',2)
hold on
plot(time_arr,y2minus,'--','Color','b','LineWidth',2)
%hold on
%patch([time_arr fliplr(time_arr)], [ y1minus fliplr(y1plus)], 'r')
hold on
plot(time_arr,y2,'-','Color','b','LineWidth',2)
hold on
box 'on'



figure(3)
y3plus = y3+hms*sqrt(y3v);
y3minus = y3-hms*sqrt(y3v);

plot(time_arr,y3_true,'s','Linewidth',2,'Color','k')
hold on
plot(time_arr,y3plus,'--','Color','g','LineWidth',2)
hold on
plot(time_arr,y3minus,'--','Color','g','LineWidth',2)
%hold on
%patch([time_arr fliplr(time_arr)], [ y1minus fliplr(y1plus)], 'r')
hold on
plot(time_arr,y3,'-','Color','g','LineWidth',2)
hold on
box 'on'



figure(4)
y4plus = y4+hms*sqrt(y4v);
y4minus = y4-hms*sqrt(y4v);

plot(time_arr,y4_true,'s','Linewidth',2,'Color','k')
hold on
plot(time_arr,y4plus,'--','Color','g','LineWidth',2)
hold on
plot(time_arr,y4minus,'--','Color','g','LineWidth',2)
%hold on
%patch([time_arr fliplr(time_arr)], [ y1minus fliplr(y1plus)], 'r')
hold on
plot(time_arr,y4,'-','Color','g','LineWidth',2)
hold on
box 'on'

figure(5)
y5plus = y5+hms*sqrt(y5v);
y5minus = y5-hms*sqrt(y5v);

plot(time_arr,y5_true,'s','Linewidth',2,'Color','k')
hold on
plot(time_arr,y5plus,'--','Color','g','LineWidth',2)
hold on
plot(time_arr,y5minus,'--','Color','g','LineWidth',2)
%hold on
%patch([time_arr fliplr(time_arr)], [ y1minus fliplr(y1plus)], 'r')
hold on
plot(time_arr,y5,'-','Color','g','LineWidth',2)
hold on
box 'on'

% figure(4)
% plot(time_arr,y4_true,'--','MarkerSize',5,'MarkerFaceColor','k')
% hold on
% plot(time_arr,y4,'-','Color','m','LineWidth',2)
% hold on
% box 'on'
% 
% figure(5)
% plot(time_arr,y5_true,'--','MarkerSize',5,'MarkerFaceColor','k')
% hold on
% plot(time_arr,y5,'-','Color','m','LineWidth',2)
% hold on
% box 'on'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y1t_list{end+1} = y1_true;
y2t_list{end+1} = y2_true;
y3t_list{end+1} = y3_true;
y4t_list{end+1} = y4_true;
y5t_list{end+1} = y5_true;
% y4t_list{end+1} = y4_true;

y1p_list{end+1} = y1;
y2p_list{end+1} = y2;
y3p_list{end+1} = y3;
y4p_list{end+1} = y4;
y5p_list{end+1} = y5;
% y4p_list{end+1} = y4;

time_list{end+1} = time_arr

end


save('validations.mat', 'y1t_list','y2t_list','y3t_list',...
     'y4t_list','y5t_list',...
     'y1p_list','y2p_list','y3p_list',...
     'y4p_list','y5p_list',...
     'y1minus','y1plus',...
     'y2minus','y2plus',...
      'y3minus','y3plus',...
      'y4minus','y4plus',...
      'y5minus','y5plus',...
     'time_list','inputs_seq_c1in', 'inputs_seq_c2in',...
     'inputs_seq_q')

