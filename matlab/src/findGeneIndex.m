function genes=findGeneIndex
load ../Data/ver3/SelGenes70.mat;
load ../Data/ver3/GeneNames.mat;        % changed since ver2
load ../Data/ver3/DrugResponse.mat;     % changed since ver2
load ../Data/ver3/GenesMutations.mat;
load ../Data/ver3/MutationCnts.mat;

GeneIndex=[];
count=0;
for i=1:1:70
    str=SelGenes70{i};
    [temp,~]=strsplit(str,'_');
    
    for j=1:1:length(temp)
        tempstr=temp{j};
        c=ismember(GeneNames,tempstr);
       % GeneIndex(i).sub(count)=find(~cellfun(@isempty,c));
       if length(find(c))
         ind=find(c);
         if(ind)
            count=count+1;
            GeneIndex.ind(count)=ind;
            GeneIndex.mut(count)=0;
            d=ismember(GenesMutations(:,1),tempstr);
            if length(find(d)) 
                 i=find(d);
                 if(i)
                    GeneIndex.mut(count)=MutationCnts(i);
                 end
            end
         end
       end
    end
end
[v,i]=sort(GeneIndex.mut,'descend');
genes=GeneIndex.ind(i);