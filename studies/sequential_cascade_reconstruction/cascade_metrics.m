function m = cascade_metrics(dsigma, centroids, true_pos, thresh)
%CASCADE_METRICS  Localization error and vertical elongation of a reconstruction.
%
%   m = cascade_metrics(dsigma, centroids, true_pos, thresh)
%
%   Quantifies the two effects the cascade is meant to improve (Section 5.8 /
%   Algorithm 1 reporting): where the reconstructed blob sits, and how
%   vertically smeared it is.  The elongation is a moment-based proxy for the
%   FWHMz/FWHMxy aspect ratio that is robust on unstructured tetrahedral
%   meshes (no gridding / interpolation required).
%
%   Inputs
%     dsigma     n x 1 reconstructed conductivity update
%     centroids  n x 3 element centroid coordinates (columns x, y, z)
%     true_pos   1 x 3 (or 3 x 1) true anomaly centre
%     thresh     support threshold as a fraction of the peak |dsigma|
%                (default 0.5); elements below it get zero weight
%
%   Output struct `m`
%     .centroid            weighted centroid of the reconstructed support
%     .localization_error  || centroid - true_pos ||
%     .rms_xy              in-plane RMS radius of the support
%     .rms_z               vertical RMS radius of the support
%     .elongation          rms_z / rms_xy  (>1 => vertically smeared)

if nargin < 4 || isempty(thresh), thresh = 0.5; end

a = abs(dsigma(:));
peak = max(a);
if peak <= eps
    % Flat image (e.g. a passing null test): metrics are undefined.
    m = struct('centroid', [NaN NaN NaN], 'localization_error', NaN, ...
               'rms_xy', NaN, 'rms_z', NaN, 'elongation', NaN);
    return;
end

wgt = a;
wgt(a < thresh * peak) = 0;
if sum(wgt) <= eps          % threshold removed everything: fall back to all
    wgt = a;
end

c = (wgt.' * centroids) / sum(wgt);        % weighted centroid (1 x 3)

d    = centroids - c;
varx = (wgt.' * (d(:,1).^2)) / sum(wgt);
vary = (wgt.' * (d(:,2).^2)) / sum(wgt);
varz = (wgt.' * (d(:,3).^2)) / sum(wgt);

m.centroid           = c;
m.localization_error = norm(c - true_pos(:).');
m.rms_xy             = sqrt(varx + vary);
m.rms_z              = sqrt(varz);
m.elongation         = m.rms_z / max(m.rms_xy, eps);
end
