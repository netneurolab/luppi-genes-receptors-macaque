%% Benchmarking macaque brain gene expression for horizontal and vertical translation.
% Authors: A.I. Luppi, Z-Q. Liu,  J.Y. Hansen,  R. Cofre,  M. Niu, E. Kuzmin,  S. Froudist-Walsh, N. Palomero-Gallagher,  & B. Misic.
 
% The study investigates the spatial correspondence of cortical patterns of gene expression
% in the macaque, against:
% (i) protein density in the macaque cortex (vertical translation); and 
% (ii) gene expression in the human cortex (horizontal translation).
 
% This repository provides code to reproduce the main results of Luppi et al., "Benchmarking macaque brain gene expression for horizontal and vertical translation." _bioRxiv_ (2024) ([preprint](https://doi.org/10.1101/2024.08.18.608440)).
 
% It was developed in MATLAB 2019a by Andrea Luppi from the the [Network Neuroscience Lab](netneurolab.github.io/) at the Montreal Neurological Institute, McGill University.
 
% This code relies on MATLAB code from the BrainSpace Toolbox(https://brainspace.readthedocs.io/en/latest/) for MATLAB by Vos de Wael et al. (2020) Communications Biology.
% The essential functions are included in this repo to ensure standalone functionality.
% For additional plotting functionality, also include in your MATLAB path the ENIGMA Toolbox
% (https://github.com/MICA-MNI/ENIGMA.git) by Lariviere et al. (2021) Nature Methods.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Repository Structure

% Main script
% The main file is Luppi_macaque_brain_gene_translation_code_4GitHub.m (this script)
% This script should work out of the box, if run from the parent directory. 
% To run, ensure you are in the main directory of the repo.

% data
% The data/ folder contains all the data you need to make this code run. 

% utils
% The utils/ folder contains support functions called by the main script, including some third-party code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET UP DIRECTORY STRUCTURE

thisDir = pwd;
addpath(genpath([thisDir, '/utils']))

materials_dir = [thisDir, '/data']
saving_dir = [thisDir, '/Results']

%Create subfolders of the saving_dir
mkdir([saving_dir, '/macaque_genes_vs_receptors/'])
mkdir([saving_dir, '/macaque_genes_vs_human_genes/'])
mkdir([saving_dir, '/gene_protein_validation/'])
mkdir([saving_dir, '/layerWise/layer_vs_layer'])
mkdir([saving_dir, '/gradients/'])

mkdir([saving_dir, '/Supplementary/stereoseq_vs_bulkRNAseq/'])
mkdir([saving_dir, '/Supplementary/macaque_genes_vs_human_RNAseq/'])
mkdir([saving_dir, '/Supplementary/cellTypes/'])
mkdir([saving_dir, '/Supplementary/structural_mediation/'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD MAIN MATERIALS

%Load macaque genes from Chen 2023 Cell, in RM atlas
load([materials_dir, '/macaque_ALL_gene_table_RM82.mat'], 'macaque_ALL_gene_table');

%Load macaque receptors from Froudist-Walsh 2023 NatureNeuro, in RM atlas
load([materials_dir, '/', 'macaque_receptors_RM82.mat'], 'receptor_table')

%Load human genes from AHBA, in human RM atlas
load([materials_dir, '/', 'human_shared_gene_table_RM82.mat'], 'human_shared_gene_table')

%Pre-computed Moran nulls from Euclidean distance of the RM atlas
load([materials_dir, '/MoranEigenvectors_RM82.mat'], 'MEM')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set some defaults
%%%%%%%%%%%%%%%%%%
%Color-schemes
brain_clrmap = 'Purples'
mat_clrmap = 'Purples' 
binary_clrmap = {[84, 40, 143]./256; [210, 210, 230]./256} %dark purple and grey

%Format for saving figures (file extension)
saving_format = 'png' %'svg' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Macaque brain genes versus receptors (ALL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Receptor name, receptor name for printing, gene IDs
receptor_gene_pairings_FULL = {
    {'AMPA'     }, {'AMPA'     }, {'GRIA1', 'GRIA2', 'GRIA3', 'GRIA4'};...
    {'kainate'  }, {'kainate'  }, {'GRIK1', 'GRIK2', 'GRIK3', 'GRIK4'};...
    {'NMDA'     }, {'NMDA'     }, {'GRIN1', 'GRIN2A', 'GRIN2B', 'GRIN2C', 'GRIN2D', 'GRIN3A', 'GRIN3B'};...
    {'GABAa'    }, {'GABA_A'    }, {'GABRA1', 'GABRA2', 'GABRA3', 'GABRA4', 'GABRA5', 'GABRA6', 'GABRB1', 'GABRB2', 'GABRB3', 'GABRD', 'GABRE', 'GABRG1', 'GABRG2', 'GABRG3', 'GABRP', 'GABRQ'};...
    {'GABAb'    }, {'GABA_B'    }, {'GABBR1', 'GABBR2'};...
    {'BZ'       }, {'GABA_{A/BZ}'       }, {'GABRA1', 'GABRA2', 'GABRA3', 'GABRA5', 'GABRB1', 'GABRB2', 'GABRB3', 'GABRD', 'GABRE', 'GABRG1', 'GABRG2', 'GABRG3', 'GABRP', 'GABRQ'};...
    {'M1'       }, {'M_1'       }, {'CHRM1'};...
    {'M2'       }, {'M_2'       }, {'CHRM2'};...
    %{'M3'       }, {'M3'       }, {'CHRM3'};... %not available 
    {'alpha1'   }, {'alpha_1'   }, {'ADRA1A', 'ADRA1B'};... ADRA1C no in macaque
    {'alpha2'   }, {'alpha_2'   }, {'ADRA2A'};... ADRA2B, ADRA2C no macaque
    {'ser_5HT1A'}, {'5HT_{1A}'}, {'HTR1A'};...
    {'ser_5HT2' }, {'5HT_{2A}' }, {'HTR2A'};... 
    {'D1'       }, {'D_1'       }, {'DRD1'};...
    };

%Set up matrix and names for the plotting function
X = [];
Y = [];
X_names = {};
Y_names = {};

for rec_num = 1:size(receptor_gene_pairings_FULL,1)
    for gen_num = 1:numel(receptor_gene_pairings_FULL{rec_num,3})
       
        x = (macaque_ALL_gene_table.(receptor_gene_pairings_FULL{rec_num,3}{gen_num}));
        y = (receptor_table.(receptor_gene_pairings_FULL{rec_num,1}{1}));

        X = [X, x];
        Y = [Y, y];

        X_names = [X_names; [receptor_gene_pairings_FULL{rec_num,3}{gen_num}]];
        Y_names = [Y_names; [receptor_gene_pairings_FULL{rec_num,2}{1}]];
        
    end
end


saving_path = [saving_dir, '/macaque_genes_vs_receptors/', 'SI_genes2receptors_ALL']

%Run correlation corrected for SA using Moran Spectral Randomisation
[corrs.macaque_genes2receptors_ALL, pvals.macaque_genes2receptors_ALL, SIG.macaque_genes2receptors_ALL,...
    zscores.macaque_genes4receptors_ALL, zscores.receptors_ALL] = ...
    fcn_correlate_and_plot_genes(X, Y, X_names, Y_names, ...
    MEM, [2,7], saving_path) 

%% FDR adjustment

%First identify which rows of the combined matrix belong to the same receptor
counter=0;
groupings_for_FDR = [];
for rec_num = 1:size(receptor_gene_pairings_FULL,1)
    counter = counter+1;
    groupings_for_FDR = [groupings_for_FDR; counter .* ones(numel(receptor_gene_pairings_FULL{rec_num,3}), 1)];
end

%Then perform FDR accounting for gene multiplicity
genes2receptors_FDR_acrossGenes_SIG = [];
genes2receptors_FDR_acrossGenes_pval = [];
for row = 1:max(groupings_for_FDR(:))

    [corrected_r, corrected_p] = fcn_fdr_matrix(corrs.macaque_genes2receptors_ALL(groupings_for_FDR==row), ...
        pvals.macaque_genes2receptors_ALL(groupings_for_FDR==row));

    genes2receptors_FDR_acrossGenes_SIG = ...
        [genes2receptors_FDR_acrossGenes_SIG; corrected_r];

    genes2receptors_FDR_acrossGenes_pval = ...
        [genes2receptors_FDR_acrossGenes_pval; corrected_p];
end

genes2receptors_FDR_acrossGenes_SIG_BIN = fcn_plot_labelled_binary_matrix(...
    genes2receptors_FDR_acrossGenes_SIG, ...
    X_names, Y_names, {''}, binary_clrmap)
saveas(gcf, [saving_dir, '/macaque_genes_vs_receptors/', 'SI_macaque_genes_vs_receptors_FDR.', 'svg'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subset used in main text of Hansen 2022 NeuroImage

receptor_gene_pairings_MAIN = {
    {'AMPA'     }, {'AMPA'     }, {'GRIA1'};...
    {'kainate'  }, {'kainate'  }, {'GRIK2'};...
    {'NMDA'     }, {'NMDA'     }, {'GRIN1'};...
    {'GABAa'    }, {'GABA_A'    }, {'GABRA1', 'GABRG2', 'GABRB2'};...
    {'GABAb'    }, {'GABA_B'    }, {'GABBR1', 'GABBR2'};...
    {'BZ'       }, {'GABA_{A/BZ}'       }, {'GABRA1', 'GABRG2', 'GABRB2'};...
    {'M1'       }, {'M_1'       }, {'CHRM1'};...
    {'M2'       }, {'M_2'       }, {'CHRM2'};...
    %{'M3'       }, {'M3'       }, {'CHRM3'};...
    {'alpha1'   }, {'alpha_1'   }, {'ADRA1A'};...
    {'alpha2'   }, {'alpha_2'   }, {'ADRA2A'};...
    {'ser_5HT1A'}, {'5HT_{1A}'}, {'HTR1A'};...
    {'ser_5HT2' }, {'5HT_{2A}' }, {'HTR2A'};...
    {'D1'       }, {'D_1'       }, {'DRD1'};...
    };

X = [];
Y = [];
X_names = {};
Y_names = {};

for rec_num = 1:size(receptor_gene_pairings_MAIN,1)
    for gen_num = 1:numel(receptor_gene_pairings_MAIN{rec_num,3})
       
        x = (macaque_ALL_gene_table.(receptor_gene_pairings_MAIN{rec_num,3}{gen_num}));
        y = (receptor_table.(receptor_gene_pairings_MAIN{rec_num,1}{1}));

        X = [X, x];
        Y = [Y, y];

        X_names = [X_names; [receptor_gene_pairings_MAIN{rec_num,3}{gen_num}]];
        Y_names = [Y_names; [receptor_gene_pairings_MAIN{rec_num,2}{1}]];
        
    end
end


saving_path = [saving_dir, '/macaque_genes_vs_receptors/', 'genes2receptors_MAIN']

%Run correlation corrected for SA using Moran Spectral Randomisation
[corrs.macaque_genes2receptors_MAIN, pvals.macaque_genes2receptors_MAIN, SIG.macaque_genes2receptors_MAIN,...
    zscores.macaque_genes4receptors_MAIN, zscores.receptors_MAIN] = ...
    fcn_correlate_and_plot_genes(X, Y, X_names, Y_names, ...
    MEM, [1,6], saving_path) 


% Store names of receptor-related genes used for main text
unique_receptor_gene_names = unique(X_names)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Brain-genes macaque vs human
% From Fulcher 2019 PNAS (plus the receptor-related genes)
%%%%%%%%%%%%%%%%%%%

%Subset of genes that are available for both human and macaque
sharedGenes_MAIN = intersect(human_shared_gene_table.Properties.VariableNames, unique_receptor_gene_names);
macaque_shared_gene_table_MAIN = macaque_ALL_gene_table(:, sharedGenes_MAIN);
human_shared_gene_table_MAIN = human_shared_gene_table(:, sharedGenes_MAIN);


%Names for plotting
for g = 1:numel(sharedGenes_MAIN)
    X_names{g,1} = [sharedGenes_MAIN{g}, ' macaque'];
    Y_names{g,1} = [sharedGenes_MAIN{g}, ' human'];
end

saving_path = [saving_dir, '/macaque_genes_vs_human_genes/', 'macaque_genes_vs_human_genes_MAIN']

%Run correlation corrected for SA using Moran Spectral Randomisation
[corrs.macaque_genes_vs_human_genes_MAIN, pvals.macaque_genes_vs_human_genes_MAIN, ...
    SIG.macaque_genes_vs_human_genes_MAIN,...
    zscores.macaque_genes_MAIN, zscores.human_genes_MAIN] = ...
    fcn_correlate_and_plot_genes(table2array(macaque_shared_gene_table_MAIN), table2array(human_shared_gene_table_MAIN), ...
    X_names, Y_names, ...
    MEM, [1,5], saving_path) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now with ALL relevant brain genes

%Subset of genes that are available for both human and macaque
sharedGenes = human_shared_gene_table.Properties.VariableNames;

%Remove those beyond Fulcher 2019 PNAS and the receptor-related ones
%(requested at revision)
load([materials_dir, '/genesNotInFulcher2019PNAS.mat'], 'genesNotInFulcher2019PNAS')
sharedGenes = setxor(sharedGenes, genesNotInFulcher2019PNAS)

%Names for plotting
for g = 1:numel(sharedGenes)
    X_names{g,1} = [sharedGenes{g}, ' macaque'];
    Y_names{g,1} = [sharedGenes{g}, ' human'];
end

saving_path = [saving_dir, '/macaque_genes_vs_human_genes/', 'SI_macaque_genes_vs_human_genes_ALL']

%Run correlation corrected for SA using Moran Spectral Randomisation
[corrs.macaque_genes_vs_human_genes_ALL, pvals.macaque_genes_vs_human_genes_ALL, ...
    SIG.macaque_genes_vs_human_genes_ALL,...
    zscores.macaque_genes_ALL, zscores.human_genes_ALL] = ...
    fcn_correlate_and_plot_genes(table2array(macaque_ALL_gene_table(:, sharedGenes)),...
    table2array(human_shared_gene_table(:, sharedGenes)), ...
    X_names, Y_names, ...
    MEM, [4,7], saving_path) 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORRELATION FOR EACH REGION - HUMAN VS MACAQUE GENES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for roi = 1:82
    [r,p] = corr(zscores.macaque_genes_ALL(roi,:)', zscores.human_genes_ALL(roi,:)', 'type', 'spearman', 'rows', 'complete');
    regional_correlation_gene2gene_rho(roi,1) = r;
    regional_correlation_gene2gene_pval(roi,1) = p;
end

%Keep track of which regions should be excluded from further analysis
roiIsNan = isnan(regional_correlation_gene2gene_rho + regional_correlation_gene2gene_pval);

save([saving_dir, '/macaque_genes_vs_human_genes/', 'Macaque_Genes_vs_Human_Genes_Regional_Corr.mat'], ...
    'regional_correlation_gene2gene_rho', 'roiIsNan')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now mRNA-protein correlation in each region
% Add other protein density maps to receptors: myelin, PV, CALB
% note that not every region will have data for every map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% T1T2 intracortical myelin marker in vivo)
load([materials_dir, '/', 'Macaque_T1T2_RM82.mat'], 'T1T2')

% Immunohistochem
load([materials_dir, '/', 'macaque_immunohisto_RM82.mat'], 'immunohisto_RM82')


%ROIS where there is no receptor data available: ignore because otherwise
%the correlation would be based on only the other 3 maps
rois2ignore = isnan(mean(zscores.receptors_ALL,2))

proteins_ALL = [zscores.receptors_ALL, ...
    T1T2, T1T2, T1T2, T1T2, ...
    immunohisto_RM82.Parvalbumin, immunohisto_RM82.Calretinin];
mRNA_ALL = [zscores.macaque_genes4receptors_ALL, ...
    macaque_ALL_gene_table.MBP, macaque_ALL_gene_table.MOBP,...
    macaque_ALL_gene_table.FTH1, macaque_ALL_gene_table.PLEKHB1,...
    macaque_ALL_gene_table.PVALB, macaque_ALL_gene_table.CALB2];

%Z-score (since not all are already z-scored)
for i = 1:size(mRNA_ALL,2)
    proteins_ALL(:,i) = (proteins_ALL(:,i)-nanmean(proteins_ALL(:,i))) ./ (nanstd(proteins_ALL(:,i)));
    mRNA_ALL(:,i) = (mRNA_ALL(:,i)-nanmean(mRNA_ALL(:,i))) ./ (nanstd(mRNA_ALL(:,i)));
end


%Correlate across gradients for each ROI 
for roi = 1:82
    [r,p] = corr(mRNA_ALL(roi,:)', proteins_ALL(roi,:)', 'type', 'spearman', 'rows', 'complete');
    regional_correlation_mRNA2protein_rho(roi,1) = r;
    regional_correlation_mRNA2protein_pval(roi,1) = p;
end
regional_correlation_mRNA2protein_rho(rois2ignore)=NaN;
regional_correlation_mRNA2protein_pval(rois2ignore)=NaN;

save([saving_dir, '/gene_protein_validation/', 'RNA_vs_protein_Regional_Corr.mat'], 'regional_correlation_mRNA2protein_rho', 'regional_correlation_mRNA2protein_pval', 'mRNA_ALL', 'proteins_ALL')

%Compare against gene-gene correlation
saving_path = [saving_dir, '/gene_protein_validation/', 'gene2gene_vs_gene2protein_regional_corrs']

%Run correlation corrected for SA using Moran Spectral Randomisation
%NB here we do not ignore negative correlations!
fcn_correlate_and_plot_genes(regional_correlation_gene2gene_rho, ...
    regional_correlation_mRNA2protein_rho, ...
    {'gene-gene correspondence'}, {'gene-protein correspondence'}, ...
    MEM, [1,1], saving_path, false)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PC1

%Find rois that are neither NaN in genes nor receptors
rois2use = find(isfinite(mean([zscores.macaque_genes4receptors_ALL, zscores.receptors_ALL], 2)));

%Brain-genes macaque
genesPC1_macaque = NaN(82,1);
[coeff,PC_score,latent] = pca(zscores.macaque_genes4receptors_ALL(rois2use, :));
genesPC1_macaque(rois2use) = PC_score(:,1);


%Brain-genes human
genePC1_human = NaN(82,1);
[coeff,PC_score,latent] = pca(zscores.human_genes_ALL(rois2use, :));
genePC1_human(rois2use) = PC_score(:,1);

%Macaque receptors PC1
receptorPC1_macaque = NaN(82,1);
[coeff,PC_score,latent] = pca(zscore(table2array(receptor_table(rois2use,:))));
receptorPC1_macaque(rois2use,1) = PC_score(:,1);


%% Compare cortical hierarchies

%Load dendritic spines from Burt2018
load([materials_dir, '/Macaque_DendriticSpines_Burt2018_RM82.mat'], 'spines_in_RM')

gradient_names_X = {'gene PC1 macaque', 'gene PC1 macaque', 'gene PC1 macaque', 'gene PC1 macaque', 'gene PC1 macaque'};
gradient_names_Y = {'T1w:T2w', 'receptor PC1 macaque', 'gene PC1 human', 'dendritic spines'}
gradients_X = [T1T2, receptorPC1_macaque, genePC1_human, spines_in_RM];
gradients_Y = [genesPC1_macaque, genesPC1_macaque, genesPC1_macaque, genesPC1_macaque];

saving_path = [saving_dir, '/gradients/', 'macaque_gradients']

%Run correlation corrected for SA using Moran Spectral Randomisation
%NB here we do not ignore negative correlations!
fcn_correlate_and_plot_genes(gradients_X, gradients_Y, gradient_names_X, gradient_names_Y, ...
    MEM, [1,4], saving_path, false)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parvalbumin and Calretinin and T1T2

% Parvalbumin and Calretinin sv PVALB and CALB2
saving_path = [saving_dir, '/gene_protein_validation/', 'genes_vs_immunohistochem']

fcn_correlate_and_plot_genes([macaque_ALL_gene_table.PVALB, macaque_ALL_gene_table.CALB2],...
    [immunohisto_RM82.Parvalbumin, immunohisto_RM82.Calretinin], ...
    {'PVALB gene macaque', 'CALB2 gene macaque'}, {'parvalbumin density', 'calretinin density'}, ...
    MEM, [1,2], saving_path)


% T1T2 vs genes from Fulcher 2019 PNAS
saving_path = [saving_dir, '/gene_protein_validation/', 'genes_vs_T1wT2w']

fcn_correlate_and_plot_genes([macaque_ALL_gene_table.MBP, macaque_ALL_gene_table.MOBP, ...
    macaque_ALL_gene_table.PVALB, macaque_ALL_gene_table.GRIN3A],...
    [T1T2, T1T2, T1T2, T1T2], ...
    {'macaque MOBP', 'macaque MBP', 'macaque PVALB', 'macaque GRIN3A'},...
    {'macaque T1w:T2w', 'macaque T1w:T2w', 'macaque T1w:T2w', 'macaque T1w:T2w'}, ...
    MEM, [1,4], saving_path, false)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Validation T1T2w from Burt
load([materials_dir, '/', 'T1T2_Burt2018_RM82.mat'], 'T1T2_Burt')

saving_path = [saving_dir, '/gene_protein_validation/', 'SI_genePC1_vs_T1T2w_Burt2018_replication']

fcn_correlate_and_plot_genes(genesPC1_macaque, T1T2_Burt, {'gene PC1 macaque'}, {'T1w:T2w (Burt 2018)'}, ...
    MEM, [1,1], saving_path, false)


%% Validation Immunohistochem from Kondo 1999

load([materials_dir, '/', 'Kondo1999_immunohistochem.mat'], 'Immunohisto_ranks_RM82')

% Parvalbumin and Calretinin sv PVALB and CALB2
saving_path = [saving_dir, '/gene_protein_validation/', 'SI_genes_vs_immunohistochem_Kondo1999_replication']

fcn_correlate_and_plot_genes([macaque_ALL_gene_table.PVALB, macaque_ALL_gene_table.CALB2],...
    [Immunohisto_ranks_RM82.PVALB, Immunohisto_ranks_RM82.CR], ...
    {'PVALB gene macaque', 'CALB2 gene macaque'}, {'parvalbumin density (rank)', 'calretinin density (rank)'}, ...
    MEM, [1,2], saving_path)





%%%%%%%%%%%%%%%%%%%%%%%
%% *LAYER-WISE* MATCH WITH RECEPTORS - ALL pairs from MURGAS 2022 NeuroImg

load([materials_dir, '/', 'macaque_GeneExpression_LayerSpecific_RM82.mat'], 'layerSpecific_GeneExpression_in_RM')

%Receptor name, receptor name for printing, gene IDs
receptor_gene_pairings_FULL = {
    {'AMPA'     }, {'AMPA'     }, {'GRIA1', 'GRIA2', 'GRIA3', 'GRIA4'};...
    {'kainate'  }, {'kainate'  }, {'GRIK1', 'GRIK2', 'GRIK3', 'GRIK4'};...
    {'NMDA'     }, {'NMDA'     }, {'GRIN1', 'GRIN2A', 'GRIN2B', 'GRIN2C', 'GRIN2D', 'GRIN3A', 'GRIN3B'};... 
    {'GABAa'    }, {'GABA_A'    }, {'GABRA1', 'GABRA2', 'GABRA3', 'GABRA4', 'GABRA5', 'GABRA6', 'GABRB1', 'GABRB2', 'GABRB3', 'GABRD', 'GABRE', 'GABRG1', 'GABRG2', 'GABRG3', 'GABRP', 'GABRQ'};...
    {'GABAb'    }, {'GABA_B'    }, {'GABBR1', 'GABBR2'};...
    {'BZ'       }, {'GABA_{A/BZ}'       }, {'GABRA1', 'GABRA2', 'GABRA3', 'GABRA5', 'GABRB1', 'GABRB2', 'GABRB3', 'GABRD', 'GABRE', 'GABRG1', 'GABRG2', 'GABRG3', 'GABRP', 'GABRQ'};...
    {'M1'       }, {'M_1'       }, {'CHRM1'};...
    {'M2'       }, {'M_2'       }, {'CHRM2'};...
    %{'M3'       }, {'M3'       }, {'CHRM3'};... %no gene from Chen 2023
    {'alpha1'   }, {'alpha_1'   }, {'ADRA1A', 'ADRA1B'};... ADRA1C no in macaque
    {'alpha2'   }, {'alpha_2'   }, {'ADRA2A'};... ADRA2B, ADRA2C no macaque
    {'ser_5HT1A'}, {'5HT_{1A}'}, {'HTR1A'};...
    {'ser_5HT2' }, {'5HT_{2A}' }, {'HTR2A'};... 
    {'D1'       }, {'D_1'       }, {'DRD1'};...
    };


for ll = 1:6

    X = [];
    Y = [];
    X_names = {};
    Y_names = {};


    for rec_num = 1:size(receptor_gene_pairings_FULL,1)
        for gen_num = 1:numel(receptor_gene_pairings_FULL{rec_num,3})

            x = (layerSpecific_GeneExpression_in_RM.(receptor_gene_pairings_FULL{rec_num,3}{gen_num}){ll});
            %zX = (x-nanmean(x)) ./ (nanstd(x));
            y = (receptor_table.(receptor_gene_pairings_FULL{rec_num,1}{1}));
            %zY = (y-nanmean(y)) ./ (nanstd(y));

            X = [X, x];
            Y = [Y, y];

            X_names = [X_names; [receptor_gene_pairings_FULL{rec_num,3}{gen_num}]];
            Y_names = [Y_names; [receptor_gene_pairings_FULL{rec_num,2}{1}]];

        end %end loop over genes
    end %end loop over receptors


    saving_path = [saving_dir, '/layerWise/', 'SI_genes2receptors_Layer', num2str(ll)]

    %Run correlation corrected for SA using Moran Spectral Randomisation
    [corrs.macaque_genes2receptors_ALL_layerWise{ll}, pvals.macaque_genes2receptors_ALL_layerWise{ll},...
        SIG.macaque_genes2receptors_ALL_layerWise{ll},...
        zscores.macaque_genes4receptors_ALL_layerWise{ll}, zscores.receptors_ALL_layerWise{ll}, ] = ...
        fcn_correlate_and_plot_genes(X, Y, X_names, Y_names, ...
        MEM, [2,7], saving_path)

end %end loop over layers


%Combined as heatmap
for ll = 1:6
    layerWise_gene_receptor_ALL_corr_mat(:,ll) = corrs.macaque_genes2receptors_ALL_layerWise{ll};
    layerWise_gene_receptor_ALL_SIG_mat(:,ll) = SIG.macaque_genes2receptors_ALL_layerWise{ll};
    layerWise_gene_receptor_ALL_pval_mat(:,ll) = pvals.macaque_genes2receptors_ALL_layerWise{ll};
end

%Summary plot with gene names on one side, receptor names on the other, layer on bottom
fcn_plot_labelled_binary_matrix(layerWise_gene_receptor_ALL_SIG_mat, X_names, Y_names, {'L1','L2', 'L3', 'L4', 'L5', 'L6' }, binary_clrmap)
saveas(gcf, [saving_dir, '/layerWise/', 'CombinedLayer_Receptors_vs_Genes_SIG_MAT.', 'svg'])


%Proportion significant
glutamate_success = mean(layerWise_gene_receptor_ALL_SIG_mat(1:15, :),1);
GABA_success = mean(layerWise_gene_receptor_ALL_SIG_mat(16:47, :),1);
neuromod_success = mean(layerWise_gene_receptor_ALL_SIG_mat(48:55, :),1); 

fcn_quick_fig([GABA_success; glutamate_success; neuromod_success], 'Proportion significant', 0, {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'}, {'GABA'; 'glutamate'; 'neuromodulators'})
try; colormap(brewermap(256, 'Purples')); end %cosmetic only
saveas(gcf, [saving_dir, '/layerWise/', 'Proportion_significant_by_ReceptorType_mat.', 'svg'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now perform FDR across genes AND layers for the same receptor

%First identify which rows of the combined matrix belong to the same receptor
counter=0;
groupings_for_FDR = [];
for rec_num = 1:size(receptor_gene_pairings_FULL,1)
    counter = counter+1;
    groupings_for_FDR = [groupings_for_FDR; counter .* ones(numel(receptor_gene_pairings_FULL{rec_num,3}), 1)];
end

%Then perform FDR accounting for both gene and layer multiplicity
layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG = [];
layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_pval = [];
for row = 1:max(groupings_for_FDR(:))

    [corrected_r, corrected_p] = fcn_fdr_matrix(layerWise_gene_receptor_ALL_corr_mat(groupings_for_FDR==row, :), ...
        layerWise_gene_receptor_ALL_pval_mat(groupings_for_FDR==row, :));

    layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG = ...
        [layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG; corrected_r];

    layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_pval = ...
        [layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_pval; corrected_p];
end

layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG_BIN = fcn_plot_labelled_binary_matrix(...
    layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG, ...
    X_names, Y_names, {'L1','L2', 'L3', 'L4', 'L5', 'L6' }, binary_clrmap)
saveas(gcf, [saving_dir, '/layerWise/', 'SI_combinedLayer_Receptors_vs_Genes_FDR_AcrossGroups_MAT.', 'svg'])

%Proportion significant
glutamate_success_FDR = mean(layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG_BIN(1:15, :),1);
GABA_success_FDR = mean(layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG_BIN(16:47, :),1);
neuromod_success_FDR = mean(layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG_BIN(48:55, :),1); %(48:57, :),1);

fcn_quick_fig([GABA_success_FDR; glutamate_success_FDR; neuromod_success_FDR], 'proportion significant after FDR', 0, {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'}, {'GABA'; 'glutamate'; 'neuromodulators'})
try; colormap(brewermap(256, 'Purples')); end % cosmetic only
saveas(gcf, [saving_dir, '/layerWise/', 'SI_proportion_significant_by_ReceptorType_FDR_mat.', 'svg'])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% stereo-seq vs bulk RNA-seq from Bo 2023 NatComm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load([materials_dir, '/', 'macaque_stereoseq_genes_sharedWithRNAseq_RM82.mat'], 'macaque_stereoseq_genes_sharedWithRNAseq')
load([materials_dir, '/', 'macaque_bulkRNAseq_genes_RM82.mat'], 'macaque_RNAseq_genes')

%Genes shared by the two modalities (since some of the genes in Bo 2023
%NatComm are not present in Chen 2023 Cell)
stereoseq_RNAseq_shared_geneIDs = macaque_stereoseq_genes_sharedWithRNAseq.Properties.VariableNames

%Prepare names for plotting
names_X = {};
names_Y = {};
for counter = 1:numel(stereoseq_RNAseq_shared_geneIDs)
    thisGene = stereoseq_RNAseq_shared_geneIDs{counter};
    names_X = [names_X; [thisGene ' stereo-seq'] ];
    names_Y = [names_Y; [thisGene ' bulk RNA-seq'] ];
end

saving_path = [saving_dir, '/Supplementary/stereoseq_vs_bulkRNAseq/', 'macaque_stereoseq_vs_bulkRNAseq']

%Run correlation corrected for SA using Moran Spectral Randomisation
[corrs.stereoseq_vs_bulk, pvals.stereoseq_vs_bulk, SIG.stereoseq_vs_bulk,...
    zscores.StereoSeq, zscores.BulkRNA] = ...
    fcn_correlate_and_plot_genes(...
    table2array(macaque_stereoseq_genes_sharedWithRNAseq(:, stereoseq_RNAseq_shared_geneIDs)), ...
    table2array(macaque_RNAseq_genes(:, stereoseq_RNAseq_shared_geneIDs)), ...
    names_X, names_Y, ...
    MEM, [1,6], saving_path) 



%% PC1
roisNoNan = not(isnan(mean([zscores.StereoSeq, zscores.BulkRNA],2)))

stereoseq_PC1_macaque = NaN(82,1);
[coeff,PC_score,latent] = pca((zscores.StereoSeq(roisNoNan, :)));
stereoseq_PC1_macaque(roisNoNan) = PC_score(:,1);

RNAseq_PC1_macaque = NaN(82,1);
[coeff,PC_score,latent] = pca((zscores.BulkRNA(roisNoNan, :)));
RNAseq_PC1_macaque(roisNoNan) = PC_score(:,1);

saving_path = [saving_dir, '/Supplementary/stereoseq_vs_bulkRNAseq/', 'stereoseq_vs_bulkRNAseq_genePC1']

X_names = {'macaque stereo-seq PC1', 'macaque receptor PC1'}
Y_names = {'macaque RNA-seq PC1', 'macaque RNA-seq PC1'}
X = [stereoseq_PC1_macaque, receptorPC1_macaque];
Y = [RNAseq_PC1_macaque, RNAseq_PC1_macaque];

fcn_correlate_and_plot_genes(X, Y, X_names, Y_names, MEM, [1,2], saving_path)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RNAseq - genes vs broadly related receptors (pertaining to same neurotransmitter)


%Receptor name, receptor name for printing, gene IDs
receptor_gene_pairings_RNAseq = {
    {'AMPA'     }, {'AMPA'     }, {'GRM4', 'GRIA4'};...
    {'kainate'  }, {'kainate'  }, {'GRM4', 'GRIA4'};...
    {'NMDA'     }, {'NMDA'     }, {'GRM4', 'GRIA4'};...
    {'GABAa'    }, {'GABA_A'    }, {'GABRR3', 'GABRQ'};...
    {'GABAb'    }, {'GABA_B'    }, {'GABRR3', 'GABRQ'};...
    {'BZ'       }, {'GABA_{A/BZ}' }, {'GABRR3', 'GABRQ'};...
    {'M1'       }, {'M_1'       }, {'CHAT', 'CHRNA1'};...
    {'M2'       }, {'M_2'       }, {'CHAT', 'CHRNA1'};...
    {'M3'       }, {'M_3'       }, {'CHAT', 'CHRNA1'};...
    {'alpha1'   }, {'alpha_1'   }, {'ADRA2C', 'PNMT'};...
    {'alpha2'   }, {'alpha_2'   }, {'ADRA2C', 'PNMT'};...
    {'ser_5HT1A'}, {'5HT_{1A}'}, {'HTR1B', 'HTR2C'};...
    {'ser_5HT2' }, {'5HT_{2A}' }, {'HTR1B', 'HTR2C'};...
    {'D1'       }, {'D_1'       }, {'DRD2', 'SLC6A2', 'NTS', 'TH'};...
    };


% Reorder to obtain the maps and names of matched genes and receptors
X = [];
Y = [];
X_names = {};
Y_names = {};

counter=0;
for rec_num = 1:size(receptor_gene_pairings_RNAseq,1)
    thisReceptor = receptor_gene_pairings_RNAseq{rec_num,1}{1};

    for gen_num = 1:numel(receptor_gene_pairings_RNAseq{rec_num,3})
        thisGene = receptor_gene_pairings_RNAseq{rec_num,3}{gen_num};

        counter = counter+1;
        x = macaque_RNAseq_genes.(thisGene);
        y = receptor_table.(thisReceptor);

        X = [X, x];
        Y = [Y, y];
        X_names = [X_names; thisGene];
        Y_names = [Y_names; receptor_gene_pairings_RNAseq{rec_num,2}{1}];

    end
end


saving_path = [saving_dir, '/Supplementary/stereoseq_vs_bulkRNAseq/', 'macaque_bulkRNAseq_vs_receptors']

%Run correlation corrected for SA using Moran Spectral Randomisation
[corrs.receptors2bulk, pvals.receptors2bulk, SIG.receptors2bulk] = ...
    fcn_correlate_and_plot_genes(X, Y, X_names, Y_names, ...
    MEM, [1,6], saving_path) 


%% FDR adjustment

%First identify which rows of the combined matrix belong to the same receptor
counter=0;
RNAseq_groupings_for_FDR = [];
for rec_num = 1:size(receptor_gene_pairings_RNAseq,1)
    counter = counter+1;
    RNAseq_groupings_for_FDR = [RNAseq_groupings_for_FDR; counter .* ones(numel(receptor_gene_pairings_RNAseq{rec_num,3}), 1)];
end

%Then perform FDR accounting for gene multiplicity
genes2receptors_RNAseq_FDR_acrossGenes_SIG = [];
genes2receptors_RNAseq_FDR_acrossGenes_pval = [];
for row = 1:max(RNAseq_groupings_for_FDR(:))

    [corrected_r, corrected_p] = fcn_fdr_matrix(corrs.receptors2bulk(RNAseq_groupings_for_FDR==row), ...
        pvals.receptors2bulk(RNAseq_groupings_for_FDR==row));

    genes2receptors_RNAseq_FDR_acrossGenes_SIG = ...
        [genes2receptors_RNAseq_FDR_acrossGenes_SIG; corrected_r];

    genes2receptors_RNAseq_FDR_acrossGenes_pval = ...
        [genes2receptors_RNAseq_FDR_acrossGenes_pval; corrected_p];
end

genes2receptors_RNAseq_FDR_acrossGenes_SIG_BIN = fcn_plot_labelled_binary_matrix(...
    genes2receptors_RNAseq_FDR_acrossGenes_SIG, ...
    X_names, Y_names, {''}, binary_clrmap)
saveas(gcf, [saving_dir, '/Supplementary/stereoseq_vs_bulkRNAseq/', 'macaque_bulkRNAseq_vs_receptors_FDR.', 'svg'])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Human microarray genes vs Macaque bulk-RNA genes

load([materials_dir, '/humanGenes_shared_with_macaqueRNAseq_RM82.mat'], 'humanGenes_shared_with_macaqueRNAseq')
macaque_bulkRNA_vs_human_shared_gene_names = humanGenes_shared_with_macaqueRNAseq.Properties.VariableNames

%Format matrices and names for the correlation function
X_names = {};
Y_names = {};

counter=0;
for i = 1: numel(macaque_bulkRNA_vs_human_shared_gene_names)
    counter = counter+1;

    thisGene = macaque_bulkRNA_vs_human_shared_gene_names{counter};

    X_names = [X_names; [thisGene ' macaque RNA-seq'] ];
    Y_names = [Y_names; [thisGene ' human'] ];

end

saving_path = [saving_dir, '/Supplementary/stereoseq_vs_bulkRNAseq/', 'macaque_bulkRNAseq_vs_human_AHBA', saving_format]

%Run correlations and save plots
[corrs.bulk2human, pvals.bulk2human, SIG.bulk2human] = ...
    fcn_correlate_and_plot_genes(table2array(macaque_RNAseq_genes(:, macaque_bulkRNA_vs_human_shared_gene_names)), ...
    table2array(humanGenes_shared_with_macaqueRNAseq(:, macaque_bulkRNA_vs_human_shared_gene_names)), X_names, Y_names, ...
    MEM, [1,4], saving_path)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Human AHBA RNA-seq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([materials_dir, '/', 'human_genes_RNAseq_in_RM82.mat'], 'human_genes_RNAseq') 

%Ensure same genes in the same order
macaque_shared_gene_table = macaque_ALL_gene_table(:, sharedGenes);
human_RNAseq_shared_gene_table = human_genes_RNAseq(:, sharedGenes);

%Names for plotting
for g = 1:numel(sharedGenes)
    X_names{g,1} = [sharedGenes{g}, ' macaque'];
    Y_names{g,1} = [sharedGenes{g}, ' human RNA-seq'];
end

saving_path = [saving_dir, '/Supplementary/macaque_genes_vs_human_RNAseq/', 'macaque_genes_vs_human_genes_RNAseq']

%Run correlation corrected for SA using Moran Spectral Randomisation
[corrs.macaque_genes_vs_human_genes_RNAseq, pvals.macaque_genes_vs_human_genes_RNAseq, ...
    SIG.macaque_genes_vs_human_genes_RNAseq,...
    ~, zscores.human_genes_RNAseq] = ...
    fcn_correlate_and_plot_genes(table2array(macaque_shared_gene_table), table2array(human_RNAseq_shared_gene_table), ...
    X_names, Y_names, ...
    MEM, [4,7], saving_path)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Human AHBA RNAseq genes vs macaque RNAseq genes

shared_genes_RNAseq = intersect(macaque_RNAseq_genes.Properties.VariableNames, human_RNAseq_shared_gene_table.Properties.VariableNames)

%Names for plotting
for g = 1:numel(shared_genes_RNAseq)
    X_names{g,1} = [shared_genes_RNAseq{g}, ' macaque RNA-seq'];
    Y_names{g,1} = [shared_genes_RNAseq{g}, ' human RNA-seq'];
end

saving_path = [saving_dir, '/Supplementary/macaque_genes_vs_human_RNAseq/', 'macaque_genes_RNAseq_vs_human_genes_RNAseq']

%Run correlation corrected for SA using Moran Spectral Randomisation
[corrs.macaque_genes_RNAseq_vs_human_RNAseq, pvals.macaque_genes_RNAseq_vs_human_RNAseq, ...
    SIG.macaque_genes_RNAseq_vs_human_RNAseq, ~, ~] = ...
    fcn_correlate_and_plot_genes(table2array(macaque_RNAseq_genes(:, shared_genes_RNAseq)), ...
    table2array(human_genes_RNAseq(:, shared_genes_RNAseq)), ...
    X_names, Y_names, ...
    MEM, [1,4], saving_path)

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Macaque gene PC1 vs human RNAseq PC1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Brain-genes human RNAseq PC1
genePC1_human_RNAseq = NaN(82,1);
[coeff,PC_score,latent] = pca(zscores.human_genes_RNAseq(rois2use, :));
genePC1_human_RNAseq(rois2use) = PC_score(:,1);

saving_path = [saving_dir, '/Supplementary/macaque_genes_vs_human_RNAseq/', 'macaque_genePC1_vs_human_RNAseq_genePC1']

%Run correlation corrected for SA using Moran Spectral Randomisation
X_names = {'macaque stereo-seq PC1', 'macaque RNA-seq PC1'}
Y_names = {'human RNA-seq PC1', 'human RNA-seq PC1'}
X = [genesPC1_macaque, RNAseq_PC1_macaque];
Y = [genePC1_human_RNAseq, genePC1_human_RNAseq];

fcn_correlate_and_plot_genes(X, Y, X_names, Y_names, MEM, [1,2], saving_path)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REGIONAL CORRELATION between receptors and RNA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for roi = 1:82
    if isnan(nanmean(zscores.human_genes_RNAseq(roi,:)))
        regional_correlation_gene2gene_RNAseq_rho(roi,1) = NaN;
    else
        [r,p] = corr(zscores.macaque_genes_ALL(roi,:)', zscores.human_genes_RNAseq(roi,:)',...
            'type', 'spearman', 'rows', 'complete');
        regional_correlation_gene2gene_RNAseq_rho(roi,1) = r;
    end
end

save([saving_dir, '/Supplementary/macaque_genes_vs_human_RNAseq/', ...
    'macaque_genes_vs_human_RNAseq_Genes_Regional_Corr.mat'], ...
    'regional_correlation_gene2gene_RNAseq_rho'); 


%% Compare inter-species correlations obtained from microarray and RNA-seq
saving_path = [saving_dir, '/Supplementary/macaque_genes_vs_human_RNAseq/', 'macaque2humanRNAseq_regional_gene_corr_']

fcn_correlate_and_plot_genes(regional_correlation_gene2gene_rho, ...
    regional_correlation_gene2gene_RNAseq_rho, ...
    {'macaque to human microarray'}, ...
   {'macaque to human RNA-seq'}, ...
   MEM, [1,1], saving_path)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CELL TYPES VS RECEPTORS

load([materials_dir, '/', 'macaque_cellTypes_norm_table_RM82.mat'], 'macaque_cellTypes_norm_table')

X = [];
Y = [];
X_names = receptor_table.Properties.VariableNames;
Y_names = macaque_cellTypes_norm_table.Properties.VariableNames;

for rc = 1: numel(X_names)
    for ct = 1:numel(Y_names)

        %Zscore manually since we don't use Moran
        x = receptor_table.(X_names{rc});
        zX = (x-nanmean(x)) ./ (nanstd(x));
        y = macaque_cellTypes_norm_table.(Y_names{ct});
        zY = (y-nanmean(y)) ./ (nanstd(y));

        %Correlate
        [cellTypes_vs_receptors_corr(rc, ct)] = ...
            corr(zX, zY, 'type', 'spearman', 'rows', 'complete');

    end
end

%Sort the matrix to highlight patterns;
%this uses BF_ClusterReorder() function from Ben Fulcher;
% if not available, loads pre-computed reorderings
try
    [sortedReceptors] = BF_ClusterReorder(cellTypes_vs_receptors_corr);
    [sortedCells] = BF_ClusterReorder(cellTypes_vs_receptors_corr(sortedReceptors, :)');
catch
    load([materials_dir, '/precomputed_cell_and_receptor_reorderings.mat'], 'sortedReceptors', 'sortedCells')
end

fcn_quick_fig(zscore(cellTypes_vs_receptors_corr(sortedReceptors, sortedCells)), ...
    ['SORTED cell types vs receptors CORR - ZSCORE'], [], Y_names(sortedCells), X_names(sortedReceptors))

try; colormap(brewermap(256, '-PRGn')); end %use brewermap to make pretty purple map if available (cosmetic only so we use try-catch to avoid breaking if unavailable)

saveas(gcf, [saving_dir, '/Supplementary/cellTypes/', 'cellTypes_vs_receptors_corr_Zscored_SORTED.', saving_format])

save([saving_dir, '/Supplementary/cellTypes/', 'cellTypes_vs_receptors_corr_Zscored_SORTED.mat'], 'cellTypes_vs_receptors_corr', 'sortedReceptors', 'sortedCells', 'X_names', 'Y_names')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORRELATION BETWEEN LAYER-SPECIFIC GENE AND RECEPTOR EXPRESSION

load([materials_dir, '/LayerSpecificGeneExpression_IPL_V1.mat'], 'layerSpecific_GeneExpression')
load([materials_dir, '/LayerSpecificReceptorDensity_IPL_V1.mat'], 'layerSpecific_ReceptorDensity')

layerSpecificRois = fieldnames(layerSpecific_GeneExpression)

for roi = 1:numel(layerSpecificRois)
    thisROI = layerSpecificRois{roi}

    X = [];
    Y = [];
    X_names = {};
    Y_names = {};

    counter=0;
    for rec_num = 1:size(receptor_gene_pairings_FULL,1)
        for gen_num = 1:numel(receptor_gene_pairings_FULL{rec_num,3})

            counter=counter+1;
            x = (layerSpecific_ReceptorDensity.(thisROI).(receptor_gene_pairings_FULL{rec_num,1}{1}));
            zX = (x-nanmean(x)) ./ (nanstd(x));
            y = (layerSpecific_GeneExpression.(thisROI).(receptor_gene_pairings_FULL{rec_num,3}{gen_num}))';
            zY = (y-nanmean(y)) ./ (nanstd(y));

            X = [X, zX];
            Y = [Y, zY];

            X_names = [X_names; [receptor_gene_pairings_FULL{rec_num,2}]];
            Y_names = [Y_names; [receptor_gene_pairings_FULL{rec_num,3}{gen_num}]];

            [gene_receptor_LayerByLayer_corr(counter,roi), gene_receptor_LayerByLayer_pval(counter,roi)] = ...
                corr(zX, zY, 'type', 'spearman');

        end %end loop over genes
    end %end loop over receptors


saving_path = [saving_dir, '/layerWise/layer_vs_layer/', 'SI_withinRegions_acrossLayers_', thisROI]

[a,b,c]=fcn_correlate_and_plot_genes(X, Y, ...
    X_names, Y_names, ...
    [], [2,7], saving_path, false)

end %end loop over rois

    
%Which gene-receptor pairs show correlated expression across layers?
glutamate_layers = 1:15;
GABA_layers = (16:47);
neuromod_layers = (48:55);

%Which gene-receptor pairs show correlated expression across layers?
glutamate_layerwise_AvgAbsCorr = mean(abs(gene_receptor_LayerByLayer_corr(glutamate_layers, :)),1);
GABA_layerwise_AvgAbsCorr = mean(abs(gene_receptor_LayerByLayer_corr(GABA_layers, :)),1);
neuromod_layerwise_AvgAbsCorr = mean(abs(gene_receptor_LayerByLayer_corr(neuromod_layers, :)),1); 

fcn_quick_fig([GABA_layerwise_AvgAbsCorr; glutamate_layerwise_AvgAbsCorr; neuromod_layerwise_AvgAbsCorr],...
    'correlation magnitude', 0,...
    layerSpecificRois, {'GABA'; 'glutamate'; 'neuromodulators'})
try; colormap(brewermap(256, 'Purples')); end %cosmetic only
saveas(gcf, [saving_dir, '/layerWise/layer_vs_layer/', 'AcrossLayers_CorrelationMagnitude_by_ReceptorType_mat.', 'svg'])

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inter-layer correlations

for roi = 1:numel(layerSpecificRois)
    thisROI = layerSpecificRois{roi}

    X = [];
    Y = [];

    counter=0;
    for rec_num = 1:size(receptor_gene_pairings_FULL,1)
        for gen_num = 1:numel(receptor_gene_pairings_FULL{rec_num,3})


            counter=counter+1;
            x = (layerSpecific_ReceptorDensity.(thisROI).(receptor_gene_pairings_FULL{rec_num,1}{1}));
            zX = (x-nanmean(x)) ./ (nanstd(x));
            y = (layerSpecific_GeneExpression.(thisROI).(receptor_gene_pairings_FULL{rec_num,3}{gen_num}))';
            zY = (y-nanmean(y)) ./ (nanstd(y));

            X = [X, zX];
            Y = [Y, zY];
        end
    end


    %Correlation of different layers against each other
    X_names = {};
    Y_names = {};
    for rl = 1:6

        X_names = [X_names; 'receptors L', num2str(rl)];
        Y_names = [Y_names; 'genes L', num2str(rl)];


        for cl = 1:6

            %Correlation of all receptors vs all genes, for pairs of layers
            [gene_receptor_LayerVSLayer_corr.(thisROI)(rl,cl),...
                gene_receptor_LayerVSLayer_pval.(thisROI)(rl,cl)] = ...
                corr(X(rl, :)', Y(cl, :)', 'type', 'spearman');

        end
    end

    fcn_quick_fig(gene_receptor_LayerVSLayer_corr.(thisROI), ...
        thisROI, [], Y_names, X_names)
    saveas(gcf, [saving_dir, '/layerWise/layer_vs_layer/', 'SI_Layer_vs_Layer_gene_receptor_correlations_ALL_', thisROI, '.', 'svg'])

    [layer_vs_layer_FDR.(thisROI), layer_vs_layer_FDR_adjustedPval.(thisROI)] = fcn_fdr_matrix(gene_receptor_LayerVSLayer_corr.(thisROI),...
        gene_receptor_LayerVSLayer_pval.(thisROI));

    fcn_quick_fig(layer_vs_layer_FDR.(thisROI), ...
        thisROI, [], Y_names, X_names)
    saveas(gcf, [saving_dir, '/layerWise/layer_vs_layer/', 'SI_Layer_vs_Layer_gene_receptor_correlations_FDR_', thisROI, '.', 'svg'])


    %Now plot the FDR-significant correlations
    for rl = 1:6
        for cl = 1:6

            if layer_vs_layer_FDR.(thisROI)(rl,cl) ~=0

                % Scatter plot
                %Heuristic to control plot size and shape
                figWidth = 200;  % Adjust base width here
                figHeight = figWidth; % Adjust height to maintain square aspect

                figure('Position', [100, 100, figWidth, figHeight]);

                %Color depends on sign
                if gene_receptor_LayerVSLayer_corr.(thisROI)(rl,cl) > 0
                    scatter(X(rl, :)', Y(cl, :)', 20, 'o', ...
                        'MarkerEdgeColor' , 'k',...
                        'MarkerFaceColor' , [220, 55, 40]./256);

                elseif gene_receptor_LayerVSLayer_corr.(thisROI)(rl,cl) < 0

                    scatter(X(rl, :)', Y(cl, :)', 20, 'o', ...
                        'MarkerEdgeColor' , 'k',...
                        'MarkerFaceColor' , [25, 105, 175]./256);
                end


                % Set axis limits to ensure data points don't lie on the axes
                xlim([min(X(rl, :)) - 0.2, max(X(rl, :)) + 0.2]);
                ylim([min(Y(cl, :)) - 0.2, max(Y(cl, :)) + 0.2]);

                % Set plot title
                if layer_vs_layer_FDR_adjustedPval.(thisROI)(rl,cl) < 0.001
                    title([' r = ', sprintf('%.2f',gene_receptor_LayerVSLayer_corr.(thisROI)(rl,cl)), '; p < 0.001'])
                elseif layer_vs_layer_FDR_adjustedPval.(thisROI)(rl,cl) < 0.01
                    title([' r = ', sprintf('%.2f',gene_receptor_LayerVSLayer_corr.(thisROI)(rl,cl)), '; p < 0.01'])
                else
                    title([' r = ', sprintf('%.2f',gene_receptor_LayerVSLayer_corr.(thisROI)(rl,cl)), ...
                        '; p = ', sprintf('%.2f',layer_vs_layer_FDR_adjustedPval.(thisROI)(rl,cl))])
                end

                xlabel(['receptors L', num2str(rl)])
                ylabel(['genes L', num2str(cl)])

                saveas(gcf, [saving_dir, '/layerWise/layer_vs_layer/', 'SI_ReceptorLayer', num2str(rl), '_vs_GeneLayer', num2str(cl), '_FDR_SIG_', thisROI, '_scatter.', 'svg'])

            end %end of whether to plot scatter

        end
    end


end %end loop over ROIs



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STRUCTURAL MEDIATION
%correlate regional receptor density with the average gene expression
%in structurally connected regions, weighted by the structural connection

%Load macaque connectome
load([materials_dir, '/macaqueTVB_82_asym.mat'])

%transpose to be in row2col format
connectome = connectome';

X = [];
Y = [];
X_names = {};
Y_names = {};

for rec_num = 1:size(receptor_gene_pairings_FULL,1)
    for gen_num = 1:numel(receptor_gene_pairings_FULL{rec_num,3})

            %Get SC-weighted average of neighbouring gene expression
            myVec = macaque_ALL_gene_table.(receptor_gene_pairings_FULL{rec_num,3}{gen_num});

            for n = 1:82
                nonZeroConnections_OUT = find(connectome(n, :) > 0);
                SC_weighted_gene_expression_OUT(n,1) = nanmean(myVec(nonZeroConnections_OUT)...
                    .* connectome(n, nonZeroConnections_OUT)' );
            end

        x = SC_weighted_gene_expression_OUT; %(geneTable.(receptor_gene_pairings{rec_num,3}{gen_num}));
        y = receptor_table.(receptor_gene_pairings_FULL{rec_num,1}{1});

        X = [X, x];
        Y = [Y, y];
        X_names = [X_names; [receptor_gene_pairings_FULL{rec_num,3}{gen_num}]];
        Y_names = [Y_names; [receptor_gene_pairings_FULL{rec_num,2}{1}]];

    end
end

saving_path = [saving_dir, '/Supplementary/structural_mediation/', 'macaque_genes_vs_receptors_SC_mediation']

fcn_correlate_and_plot_genes(X, Y, X_names, Y_names, ...
   MEM, [2,7], saving_path)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPARE PROPORTIONS

glutamate_success_vec = layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG_BIN(glutamate_layers, :);
glutamate_success_vec = glutamate_success_vec(:);

GABA_success_vec = layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG_BIN(GABA_layers, :);
GABA_success_vec = GABA_success_vec(:);

neuromod_success_vec = layerWise_gene_receptor_ALL_FDR_acrossLayersAndGenes_SIG_BIN(neuromod_layers, :);
neuromod_success_vec = neuromod_success_vec(:);

%Chi2
[chi2_stats_layerWise.glut_GABA.chi2, ...
    chi2_stats_layerWise.glut_GABA.p, ...
    chi2_stats_layerWise.glut_GABA.summary] = ...
    fcn_chi2_test(glutamate_success_vec, GABA_success_vec, false)

[chi2_stats_layerWise.glut_neuromod.chi2, ...
    chi2_stats_layerWise.glut_neuromod.p, ...
    chi2_stats_layerWise.glut_neuromod.summary] = ...
    fcn_chi2_test(glutamate_success_vec, neuromod_success_vec, false)

[chi2_stats_layerWise.GABA_neuromod.chi2,...
    chi2_stats_layerWise.GABA_neuromod.p, ...
    chi2_stats_layerWise.GABA_neuromod.summary] = ...
    fcn_chi2_test(GABA_success_vec, neuromod_success_vec, false)
save([saving_dir, '/layerWise/', 'chi2_stats_layerWise.mat'], 'chi2_stats_layerWise')
