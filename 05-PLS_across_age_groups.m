
% PLS code (pls_analysis.m) can be downloaded at
% http://pls.rotman-baycrest.on.ca/source/ ("Latest PLS Applications")

%% load
clear all
clc
cd '~/Documents/SickKids/abagen_analysis/'
load('GeneExpression.mat')           % relevant node indices
load('SpinsTwirls.mat')           % spin test indices
load('coordinates.mat')          % (x,y,z) coordinates for brain regions
load('GeneNames.mat')          % (x,y,z) coordinates for brain regions

corridinates= table2array(coordinates);
GeneExpression= table2array(GeneExpressionDestrieux);

SpinsTwrils=SpinsTwrils+1;

addpath('./natsortfiles')

topgene=readtable('~/Documents/SickKids/abagen_analysis/OverlapGeneLoadings.csv');

%% set up of PLS analysis

addpath(genpath('./Pls/'));

% set up PLS analysis

X = zscore(GeneExpression);

nnodes = 148; % number of nodes/ ROIs 
ngenes = length(GeneExpression);
nterms= 4;
bootnnodes = 118; % bootstrapped number of ROIs

% behav pls
option.method = 3;
option.num_boot = 0;
option.num_perm = 0;               % zero permutations because they will be run manually later to account for spatial autocorrelation

exp{1} = X;

filename= dir('/Users/jason/Documents/SickKids/ICCvalues2/ICC_group_mean_*_arrhythmic.csv');

sort_files=natsort({filename.name})';

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

%% iterate over age groups

for i=1:37
    
    filename2PLS= [filename(i).folder, '/', sort_files{i}];
    
    ICC=readmatrix(filename2PLS);
    X = zscore(GeneExpression);
    exp{1} = X;

    Y = zscore(ICC);
    option.stacked_behavdata = Y;

    result = pls_analysis(exp, nnodes, 1, option); % this is the PLS result that is used in all other analyses

    cb=(result.s.^2)/(sum(result.s.^2)); % calculate percent varaince explained 

    var_explained(i,:)=cb;
    
    %% spin test
    % this code comes from pls_analysis.m and is modified to account for a
    % spatial autocorrelation-preserving permutation test

    Y = zscore(ICC);
    nspins = 1000;                 % number of permutations ("spins")
    s_spins = zeros(nterms,nspins); % singular values
    option.method = 3;              % set up PLS
    option.num_boot = 0;
    option.num_perm = 0;
    X = zscore(GeneExpression);
    exp{1} = X;
    for k = (1+ (i-1)*1000):(nspins*i)    
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

    pvalues(i,:)=sprob;
    
    var_CI=[];
    for b=1:1000

        bootind= randi(148,[1,bootnnodes]);
        % set up PLS analysis

        X = zscore(GeneExpression(bootind,:));
        Y = zscore(ICC(bootind,:));
        option.stacked_behavdata = Y;

        exp{1} = X;

        result_CI = pls_analysis(exp, bootnnodes, 1, option); % this is the PLS result that is used in all other analyses
        cb_CI=(result_CI.s.^2)/(sum(result_CI.s.^2)); % calculate percent varaince explained 
        var_CI(b)=cb_CI(1);

    end

    var_explained_CI(i,:)=quantile(var_CI, [0.025, 0.975]);
    
    % get gene sets
    % compute the loading of each gene as the correlation between the original
    % data and the gene scores
    gload = zeros(ngenes,1);
    for k = 1:ngenes
        gload(k) = corr(GeneExpression(:,k),result.vsc(:,1));
    end
   
    overlap(i)=abs(corr(topgene.loadings, gload(topgene.index)));
    
    if corr(topgene.loadings, gload(topgene.index)) > 0
        neurophyscores(i,:)=-1*corr(ICC,result.usc(:,1)); 
    else corr(topgene.loadings, gload(topgene.index)) < 0
        neurophyscores(i,:)=corr(ICC,result.usc(:,1)); 
    end
    

end

pvalFDR2= mafdr(pvalues(:,1),'BHFDR',true);

varNames = {'varexplained', 'CIlower','CIupper','pval','pvalFDR','overlapGeneLoad', 'NeuroPhysTheta', 'NeuroPhysAlpha', 'NeuroPhysBeta', 'NeuroPhysGamma'};
Tableoutput=table(var_explained(:,1), var_explained_CI(:,1),var_explained_CI(:,2),pvalues(:,1), pvalFDR2, overlap', neurophyscores(:,1),neurophyscores(:,2),neurophyscores(:,3),neurophyscores(:,4), 'VariableNames',varNames);

writetable(Tableoutput, './PLS_genes_across_lifespan_2.csv');
