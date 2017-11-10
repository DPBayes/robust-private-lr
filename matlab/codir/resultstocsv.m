% Combine and convert the result data as csv files
% Change filenames (corr/wpc,e1/e2) if needed!
% cd where needed
load('result.mat');

% PLR: private linear regression
a = result.plr;
sn = size(a,1);
cv = size(a,2);
res = NaN(cv,sn);
for i = 1:sn
    res(:,i) = [a{i,:}];
end
m = nanmean(res);
s = std(res);
csvwrite('corr-e2-plr-mean.csv',m);
csvwrite('corr-e2-plr-std.csv',s);

% OPLR: output perturbed linear regression
a = result.oplr;
sn = size(a,1);
cv = size(a,2);
res = NaN(cv,sn);
for i = 1:sn
    res(:,i) = [a{i,:}];
end
m = nanmean(res);
s = std(res);
csvwrite('corr-e2-oplr-mean.csv',m);
csvwrite('corr-e2-oplr-std.csv',s);

% FMLR: functional mechanism linear regression
a = result.fmlr;
sn = size(a,1);
cv = size(a,2);
res = NaN(cv,sn);
for i = 1:sn
    res(:,i) = [a{i,:}];
end
m = nanmean(res);
s = std(res);
csvwrite('corr-e2-fmlr-mean.csv',m);
csvwrite('corr-e2-fmlr-std.csv',s);
