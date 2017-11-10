%path to the CCLE dataset alternative to Sanger
path='/users/mrinal/projects/pm/DATA/broad/CTRP_DR_GE_Data_Processed_20160415/missing_cases';

%load overlapping information
load([path,'/UsefulGenesIndexes.mat']);
load([path,'/MatchingDrugIndexes.mat']);

%load data
load([path,'/GeneExpression.mat']);
load([path,'/DrugResponse.mat']);
%{
load([path,'/UsefulGenes.mat']);
ccle_genes=UsefulGenes;
load(['ver2/UsefulGenes.mat']);
sanger_genes=UsefulGenes;
%}

% take the overlap
X=GeneExpression(:,UsefulGenesIndexes); %n x d
Y=DrugResponse(:,sangerDrugIndex); % n x m

%run lasso
[n,d]=size(X);
[n,m]=size(Y);

m=1;
for j=1:m
    y=Y(:,j);
    b=lasso(X,y); %,'Lambda',1);
    B(:,j)=mean(b,2);
end
avg_b=mean(B,2);
[val,ind]=sort(avg_b,'descend');


load(['ver2/UsefulGenesIndexes.mat']);
sorted_genes=UsefulGenesIndexes(ind);

save(['ver2/lasso_sorted_genes.mat'],'sorted_genes');


    

