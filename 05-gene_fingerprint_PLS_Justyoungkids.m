
% PLS code (pls_analysis.m) can be downloaded at
% http://pls.rotman-baycrest.on.ca/source/ ("Latest PLS Applications")

%% load
clear all
clc
cd '~/Documents/SickKids/abagen_analysis/'
load('ICCUnder8.mat')      % node by ICC
load('GeneExpression.mat')           % relevant node indices
load('Spins.mat')           % spin test indices
load('coordinates.mat')          % (x,y,z) coordinates for brain regions
load('GeneNames.mat')          % (x,y,z) coordinates for brain regions

corridinates= table2array(coordinates);
ICC= table2array(ICCunder8);
GeneExpression= table2array(GeneExpressionDestrieux);

clear GeneExpressionDestrieux ICCunder8

%% PLS analysis

addpath(genpath('./Pls/'));

% set up PLS analysis

X = zscore(GeneExpression);
Y = zscore(ICC);

nnodes = 148; % number of nodes/ ROIs 
ngenes = length(GeneExpression);
nterms= 4;

% behav pls
option.method = 3;
option.num_boot = 1000;
option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
option.stacked_behavdata = Y;

exp{1} = X;

result = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses
save('./result_under8.mat','result', '-v7.3')


load('./result_under8.mat')


%% spin test
% this code comes from pls_analysis.m and is modified to account for a
% spatial autocorrelation-preserving permutation test

nspins = 1000;                 % number of permutations ("spins")
s_spins = zeros(nterms,nspins); % singular values
option.method = 3;              % set up PLS
option.num_boot = 0;
option.num_perm = 0;
exp{1} = X;
for k = 1:nspins    
    option.stacked_behavdata = Y(SpinsTwrils(:,k),:);  % permute neurosynth matrix
    
    datamatsvd=rri_xcor(option.stacked_behavdata,exp{1},0); % refer to pls_analysis.m
    [r,c] = size(datamatsvd);
    if r <= c
        [pu, sperm, pv] = svd(datamatsvd',0);
    else
        [pv, sperm, pu] = svd(datamatsvd,0);
    end
    
    %  rotate pv to align with the original v
    rotatemat = rri_bootprocrust(result.v,pv);
 
    %  rescale the vectors
    pv = pv * sperm * rotatemat;

    sperm = sqrt(sum(pv.^2));
    
    s_spins(:,k) = sperm;
end

sprob = zeros(nterms,1); % p-value for each latent variable

for k = 1:nterms % get permuted (via spin test) p-values
    sprob(k) = (1+(nnz(find(s_spins(k,:)>=result.s(k)))))/(1+nspins);
end  

% just on the threshold of significance 9.9900e-04


%% plot variance explained

cb=(result.s.^2)/(sum(result.s.^2)); % calculate percent varaince explained 
cb_spin=(s_spins.^2)./repmat((sum(s_spins.^2, 1)),[4,1]); % calculate percent varaince explained 


figure
hold on
boxplot(cb_spin'*100)
plot(1:4,cb*100,'b.','MarkerSize', 30); % plot percent var explained
hold on
set(findall(gcf,'-property','FontSize'),'FontSize',12)
plot(1:4,cb*100,'b-','LineWidth', 1.5); % plot percent var explained
xlabel("Component Number")
ylabel("Percent Covariance Explained (%)")
ylim([0 100])


saveas(gcf,'./NeuroPhys_percVarEx.png')
saveas(gcf,'./NeuroPhys_percVarEx.pdf')
%saveas(gcf,'./NeuroPhys_percVarEx.fig')

% first significant component explaines 67.39% of variance 


%% bootstrap the data to set CI for the var explained 

nnodes = 118; % number of nodes/ ROIs 

var_explained_CI=[];
for b=1:1000
    
    bootind= randi(148,[1,nnodes]);
    % set up PLS analysis

    X = zscore(GeneExpression(bootind,:));
    Y = zscore(ICC(bootind,:));
    % behav pls
    option.method = 3;
    option.num_boot = 0;
    option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation
    option.stacked_behavdata = Y;

    exp{1} = X;

    result_CI = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses
    cb_CI=(result_CI.s.^2)/(sum(result_CI.s.^2)); % calculate percent varaince explained 
    var_explained_CI(b)=cb_CI(1);

end

quantile(var_explained_CI, [0.025, 0.975])

% ver explained CI 40.93% and 80.4%



%% loadings of bands

l1=corr(ICC,result.usc(:,1)); 

bands= categorical({'2.theta', '3.alpha', '4.beta', '5.gamma'});
bands = reordercats(bands,cellstr(bands)');

figure
clear g
g(1,1)=gramm('x',bands,'y',l1, 'color', bands);
g(1,1).stat_summary('geom','bar','setylim',true);
g(1,1).set_title('Neurophysiology loadings ''geom'',''bar''');

g.draw();

saveas(gcf,'./NeuroPhys_loadings.png')
saveas(gcf,'./NeuroPhys_loadings.pdf')
saveas(gcf,'./NeuroPhys_loadings.fig')

%% get genes with entrezID

T = table2cell(readtable('gene_entrez_ids')); % load entrezID of genes

gene_name = GeneNames;     % get relevant gene names
entrezIDs = zeros(size(gene_name));

idx = [];
for k = 1:length(gene_name)                                                % for each gene
    if ismember(gene_name{k}, T(:,1))                                      % if the gene has an entrezID
        entrezIDs(k) = cell2mat(T(find(strcmp(gene_name{k}, T(:,1))),2));  % store the entrezID
        idx = [idx;k];                                                     % also store the index of the gene
    end
end
%entrezIDs = entrezIDs(entrezIDs ~= 0);                                     % remove all genes without entrezID
entrezIDsNONID = entrezIDs(entrezIDs ~= 0); % this will be our background genes to compare to in the enrichment analysis 

%% get category scores

% get gene sets
% compute the loading of each gene as the correlation between the original
% data and the gene scores
gload = zeros(ngenes,1);
for k = 1:ngenes
    gload(k) = corr(GeneExpression(:,k),result.vsc(:,1));
end

ipos = find(gload > 0); % index of genes with positive loading
ineg = find(gload < 0); % index of genes with negative loading
gload_pos = gload(gload > 0); % loading of genes with positive loading
gload_neg = gload(gload < 0); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos); % sorted
[~,Ineg] = sort(gload_neg); % sorted 

threshold = 0.5; % top 50% of pos/neg genes constitute each gene set

gpos_idx = Ipos(end-floor(threshold*length(gload_pos)):end); % top 50% of genes with positive loading
gneg_idx = Ineg(1:floor(threshold*length(gload_neg)));       % top 50% of genes with negative loading

gpos_ID = entrezIDs(ipos(gpos_idx));   % these are the entrezIDs of the genes in the positive set
gneg_ID = entrezIDs(ineg(gneg_idx));  % these are the entrezIDs of the genes in the negative set

gpos_ID_nonzero= gpos_ID(gpos_ID~=0);
gneg_ID_nonzero= gneg_ID(gneg_ID~=0);

% find genes with no ids
% compute category score for them too 
ids_pos=ipos(gpos_idx);
ids_neg=ineg(gneg_idx);


% %% make a nice table with all info 
% 
varNames = {'gene', 'EntrezID','loading'};
NEGATIVE_tab=table(gene_name(ids_neg), entrezIDs(ids_neg), gload(ids_neg), 'VariableNames',varNames);

POSITIVE_tab=table(gene_name(ids_pos),entrezIDs(ids_pos), gload(ids_pos), 'VariableNames',varNames);

writetable(POSITIVE_tab, './Positive_gene_loadings.csv');
writetable(NEGATIVE_tab, './Negative_gene_loadings.csv');



%% Gene Ontology Plot

enrichmentallGOpositive=importfile('~//Documents/SickKids/abagen_analysis/SickKidsYoungCohort_positive_loadings.csv');
enrichmentallGOpositiveClean=enrichmentallGOpositive(enrichmentallGOpositive.EnrichmentFDR <0.05, :);

Loadings=zeros(length(enrichmentallGOpositiveClean.Genes),1);
for i =1:length(enrichmentallGOpositiveClean.Genes)
    
    genes_temp=split(enrichmentallGOpositiveClean.Genes(i));
    temp_gload= nan(length(genes_temp),1);
    for k =1:length(genes_temp)

        if ~isempty(find(strcmp(genes_temp(k), string(GeneNames))))
        
            temp_gload(k,1)=gload(find(strcmp(genes_temp(k), string(GeneNames))));
        end
    end
   Loadings(i)= mean(temp_gload, 'omitnan' );
end

enrichmentallGOpositiveClean.Loadings=Loadings;

cutoff=sort(enrichmentallGOpositiveClean.EnrichmentFDR, "ascend");
enrichmentallGOpositiveClean2= enrichmentallGOpositiveClean(enrichmentallGOpositiveClean.EnrichmentFDR <=cutoff(35), : );
[B ,I]= sort(enrichmentallGOpositiveClean2.FoldEnrichment, 'descend');
enrichmentallGOpositiveClean2= enrichmentallGOpositiveClean2(I,:);


m1 = floor(70*0.5);
r = (0:m1-1)'/max(m1,1);
g = r;
r = [r; ones(m1+1,1)];
g = [g; 1; flipud(g)];
b = flipud(r);
c = [r g b]; 

enrichmentallGOpositiveClean2.Pathway=lower(enrichmentallGOpositiveClean2.Pathway);


figure
wordcloud(enrichmentallGOpositiveClean2,'Pathway','Loadings'); %, 'Color', c(71:-1:37,:));
title("gene process (positive)")

saveas(gcf,'./wordclous_positive_GO.png')
saveas(gcf,'./wordclous_positive_GO.pdf')
saveas(gcf,'./wordclous_positive_GO.fig')


% negative 
enrichmentallGOnegative=importfile('~/Documents/SickKids/abagen_analysis/SickKidsYoungCohort_negative_loadings.csv');

enrichmentallGOnegativeClean=enrichmentallGOnegative(enrichmentallGOnegative.EnrichmentFDR <0.05, :);

Loadings=zeros(length(enrichmentallGOnegativeClean.Genes),1);
for i =1:length(enrichmentallGOnegativeClean.Genes)
    
    genes_temp=split(enrichmentallGOnegativeClean.Genes(i));
    temp_gload= nan(length(genes_temp),1);
    for k =1:length(genes_temp)

        if ~isempty(find(strcmp(genes_temp(k), string(GeneNames))))
        
            temp_gload(k,1)=gload(find(strcmp(genes_temp(k), string(GeneNames))));
        end
    end
   Loadings(i)= mean(temp_gload, 'omitnan' );
end


enrichmentallGOnegativeClean.Loadings=Loadings;

cutoff=sort(enrichmentallGOnegativeClean.EnrichmentFDR, "ascend");
enrichmentallGOnegativeClean2= enrichmentallGOnegativeClean(enrichmentallGOnegativeClean.EnrichmentFDR <=cutoff(35), : );
[B ,I]= sort(enrichmentallGOnegativeClean2.FoldEnrichment, 'descend');
enrichmentallGOnegativeClean2= enrichmentallGOnegativeClean2(I,:);
enrichmentallGOnegativeClean2.Loadings = abs(enrichmentallGOnegativeClean2.Loadings);



m1 = 70*0.5;
r = (0:m1-1)'/max(m1-1,1);
g = r;
r = [r; ones(m1,1)];
g = [g; flipud(g)];
b = flipud(r);
c = [r g b]; 

enrichmentallGOnegativeClean2.Pathway=lower(enrichmentallGOnegativeClean2.Pathway);

figure
wordcloud(enrichmentallGOnegativeClean2,'Pathway','FoldEnrichment' );%, 'Color', c(35:-1:1,:));
title("gene process (negative)")

saveas(gcf,'./wordclous_negative_GO.png')
saveas(gcf,'./wordclous_negative_GO.pdf')
saveas(gcf,'./wordclous_negative_GO.fig')



