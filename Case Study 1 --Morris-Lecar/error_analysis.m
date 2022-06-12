clear all
clearvars

nreps = 5;
ndesign = 12;
nnoise=1;
%dataps = [30,60,90,120,150,180];
dataps = [20,40,60,80,100,120,140,160,180,200,220,240]
y1ve = zeros(nreps, ndesign,nnoise);
y2ve = zeros(nreps, ndesign,nnoise);
mte = zeros(nreps, ndesign,nnoise);

noise =[0.0,0.0];

for l = 1:nnoise
for j=1:ndesign
    for i=1:nreps 
        mte(i,j,l) = get_mte(dataps(j),noise(l,:)) ;
    end
end
end

save('errors.mat','mte')