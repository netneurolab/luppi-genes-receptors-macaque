function[chi2stat, p, word_summary, df] = fcn_chi2_test(vec1, vec2, flag_smallSampleContinuityCorrection)
% Chi-square test to compare two proportions
% (wrapper around the MATLAB function prop_test)
% INPUTS:
% vec1 and vec2 must be binary vectors of successes (1) and failures (0);
% they can have different numbers of entries
% flag_smallSampleContinuityCorrection (boolean, default is false): 
% use Yates continuity correction for small samples?
%
% OUTPUTS:
% chi2stat (scalar): the chi-square test statistic
% p (scalar): p-value
% word_summary: a text string to facilitate reporting of the results,
% including the proportions of successes in each vector
% df (scalar): the test's degrees of freedom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smallSampleContinuityCorrection is true by default
if nargin == 2
    flag_smallSampleContinuityCorrection = false;
end

% Perform the chi-squared test
[h, p, chi2stat,df] = prop_test([sum(vec1), sum(vec2)], [numel(vec1),...
    numel(vec2)], flag_smallSampleContinuityCorrection); % chi2gof(contingencyTable, 'ctable');

word_summary = ['Vec1: ', num2str(sum(vec1)), ' out of ', num2str(numel(vec1)),...
    ' (', num2str(100 * sum(vec1)/numel(vec1)), '%); ', ...
    'Vec2: ', num2str(sum(vec2)), ' out of ', num2str(numel(vec2)),...
    ' (', num2str(100.*sum(vec2)/numel(vec2)), '%); ', ...
    'Chi-squared statistic = ', sprintf('%.2f',chi2stat),...
    '; p-value = ', sprintf('%.3f',p)]
