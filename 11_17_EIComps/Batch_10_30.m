% Batch 10.30

[ALLTS,meansC,stC,stMeans,stMeds] = scoreIFs;
save('SWData_10_30.mat');
clear
[ALLTS2,meansC2,stC2,stMeans2,stMeds2] = scoreIFs2;
save('LWData_10_30.mat');