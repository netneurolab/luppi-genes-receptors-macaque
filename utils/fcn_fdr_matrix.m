function[fdr_significant_results, fdr_pvals] = fcn_fdr_matrix(raw_results_mat, raw_pvals_mat, alpha);
% Perform filtering of results by False Discovery Rate correction
% significance
%
% INPUTS:
% raw_results_mat: R-by-C matrix of test statistics (for example, the
% correlation coefficient rho)
% raw_pvals_mat: R-by-C matrix of uncorrected p-values
% alpha (scalar): the alpha level for significance (optional: default 0.05)
%
% OUTPUTS:
% fdr_significant_results: R-by-C matrix of the results, where any entries
% that do not pass significnce after FDR correction are set to zero
% fdr_pvals:  R-by-C matrix of the FDR-corrected adjusted p-values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end

if numel(raw_results_mat) == 1

    %Deal with the case where there is only one entry (no correction)
    fdr_pvals = raw_pvals_mat;
    fdr_significant_results = mask_matrix(raw_results_mat, fdr_pvals < alpha);

else

    if (issymmetric(raw_results_mat)) % if the matrix is symmetrix, only use the upper triangular

        vec = helper_triu2vec(raw_pvals_mat);

        % Apply FDR correction with Benjamini-Hochberg method
        [fdr_mask_vec, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(vec, alpha); 
        fdr_mask_matrix = helper_vec2mat(fdr_mask_vec, size(raw_results_mat,1));


    else
        % Apply FDR correction with Benjamini-Hochberg method
        [fdr_mask_matrix, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(raw_pvals_mat, alpha); 

    end

    fdr_significant_results = mask_matrix(raw_results_mat, fdr_mask_matrix);
    fdr_pvals = adj_p; 

end

%% HELPER FUNCTIONS

% Helper function to mask a matrix A with a second matrix B
function [masked_matrix] = mask_matrix(original_matrix, mask)
    masked_matrix = zeros(size(original_matrix));
    include = find(mask==1);
    masked_matrix(include) = original_matrix(include);
end

%Takes a matrix and turns its triu (excluding diagonal) into a vector
function[vec] = helper_triu2vec(mat)
    vec = mat(not(tril(true(size(mat)))));
end

end %EOF