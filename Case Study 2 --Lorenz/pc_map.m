function [pce,elapsedtime] = pc_map(X,Y,Ns,Np,Nd)
clc;

minvals = min(X,[],1);
maxvals = max(X,[],1);
fact = 1.5;
lf = 0.95;
uf = 1.05;

x1m = mean(X(:,1));
x1s = std(X(:,1));

x2m = mean(X(:,2));
x2s = std(X(:,2));

x3m = mean(X(:,3));
x3s = std(X(:,3));


inputOpts.Marginals(1).Type = 'Gaussian';
inputOpts.Marginals(1).Parameters= [x1m x1s];

inputOpts.Marginals(2).Type = 'Gaussian';
inputOpts.Marginals(2).Parameters = [x2m x2s];

inputOpts.Marginals(3).Type = 'Gaussian';
inputOpts.Marginals(3).Parameters = [x3m x3s];



for jj=(Ns+1):(Ns+Np)
inputOpts.Marginals(jj).Type = 'Uniform';
inputOpts.Marginals(jj).Parameters= [lf*min(X(:,jj)),uf*max(X(:,jj))];
end


myInput = uq_createInput(inputOpts);
Ntrain  = Nd;
Ndata = size(X,1);
Nvp = 50;

XT = X(1:Ntrain,:);
XV = X((Ntrain+1):(Ntrain+Nvp),:);

YT = Y(1:Ntrain,:);
YV = Y((Ntrain+1):(Ntrain+Nvp),:);

%load('input_uniform.mat')
%full_input = uq_createInput(Uinput);
Nout = size(Y,2) -0;%How many outputs in the (reduced) target space 
%% Build the surrogate for the continuous responses  
fprintf('Beggining PCE calibration \n');
lab_ind = 0; %Initialize index for label
for i= 1:Nout
tic
fprintf('-------------------------------------- \n');
fprintf('Estimating for label %d .... .\n ', i);
fprintf('-------------------------------------- \n');
%lab_ind = lab_ind + 1; %Loops over labels

CMetaOpts = []; %Initialize MetaOpts field
CMetaOpts.Type = 'Metamodel';
CMetaOpts.MetaType = 'PCK';
CMetaOpts.Method = 'LARS';
CMetaOpts.Degree = 1:5;
%CMetaOpts.Corr.Family = 'Linear'
% CMetaOpts.LARS.LarsEarlyStop = false;
CMetaOpts.Display = 1;
CMetaOpts.TruncOptions.qNorm = 0.7:0.05:0.85;
% CMetaOpts.TruncOptions.MaxInteraction = 15 ;
CMetaOpts.ExpDesign.X = XT;
CMetaOpts.ExpDesign.Y = YT(:,i);
CMetaOpts.ValidationSet.X = XV;
CMetaOpts.ValidationSet.Y = YV(:,i);
pce{i,1}= uq_createModel(CMetaOpts);
elapsedtime(i)=toc;
end


end