function [p, h, stats] = mysignrank(x,y,alpha)
%SIGNRANK Wilcoxon signed rank test for zero median.
%   P = SIGNRANK(X) performs a two-sided signed rank test of the hypothesis
%   that the data in the vector X come from a distribution whose median
%   (and mean, if it exists) is zero, and returns the p-value from the
%   test.  P is the probability of observing the given result, or one more
%   extreme, by chance if the null hypothesis ("median is zero") is true.
%   Small values of P cast doubt on the validity of the null hypothesis.
%   The data are assumed to come from a continuous distribution, symmetric
%   about its median.
%
%   P = SIGNRANK(X,M) performs a two-sided test of the hypothesis that the
%   data in the vector X come from a distribution whose median is M.  M
%   must be a scalar.
%
%   P = SIGNRANK(X,Y) performs a paired, two-sided test of the hypothesis
%   that the difference between the matched samples in the vectors X and Y
%   comes from a distribution whose median is zero.  The differences X-Y
%   are assumed to come from a continuous distribution, symmetric about its
%   median.  X and Y must be the same length.
%
%   [P,H] = SIGNRANK(...) returns the result of the hypothesis test,
%   performed at the 0.05 significance level, in H.  H==0 indicates that
%   the null hypothesis ("median is zero") cannot be rejected at the 5%
%   level. H==1 indicates that the null hypothesis can be rejected at the
%   5% level.
%
%   [P,H] = SIGNRANK(...,ALPHA) returns the result of the hypothesis test
%   performed at the significance level ALPHA.
%
%   [P,H,STATS] = SIGNRANK(...) returns STATS, a structure with one or two
%   fields.  The field 'signedrank' contains the value of the signed rank
%   statistic.  If the sample size is large, then P is calculated using a
%   normal approximation and the field 'zval' contains the value of the
%   normal (Z) statistic.
%
%   See also SIGNTEST, RANKSUM, TTEST, ZTEST.

%   References:
%      [1] Hollander, M. and D. A. Wolfe.  Nonparametric Statistical
%          Methods. Wiley, 1973.
%      [2] Gibbons, J.D.  Nonparametric Statistical Inference,
%          2nd ed.  M. Dekker, 1985.

%   Copyright 1993-2003 The MathWorks, Inc. 
%   $Revision: 1.16 $

alpha=0.05;

if (length(alpha)>1)
   error('SIGNRANK requires a scalar ALPHA value.');
end
if (isnan(alpha) | (alpha <= 0) | (alpha >= 1))
   error('SIGNRANK requires 0 < ALPHA < 1.');
end

[rowx, colx] = size(x);
if (length(y) == 1)
   y = repmat(y, rowx, colx);
end
[rowy, coly] = size(y);

if min(rowx, colx) > 1 | min(rowy,coly) > 1,
   error('SIGNRANK requires vector rather than matrix data.');
end 
if rowx == 1
   rowx = colx;
   x = x';
end
if rowy == 1,
   rowy = coly;
   y = y';
end
   
if rowx ~= rowy,
   error('SIGNRANK requires the data vectors to have the same number of elements.');
end

diffxy = x - y;

diffxy(isnan(diffxy)) = [];
if (length(diffxy)==0), error('No data remaining after removal of NaNs.'), end

nodiff = find(diffxy == 0);
diffxy(nodiff) = [];
n = length(diffxy);

if (n == 0)         % degenerate case, all ties
   p = 1;
   if (nargout > 1)
      h = 0;
      if (nargout > 2)
         stats.signedrank = 0;
      end
   end
   return
end

% Find negative differences and ranks of absolute differences
neg = find(diffxy<0);
[tierank, tieadj] = tiedrank(abs(diffxy));

% Compute signed rank statistic (most extreme version)
w = sum(tierank(neg));
w = min(w, n*(n+1)/2-w);

% Compute significance exactly or with a normal approximation
if n > 15,
   z = (w-n*(n+1)/4) / sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
   p = 2*normcdf(z,0,1);
   if (nargout > 2)
      stats.zval = z;
   end
else
   allposs = (ff2n(n))';
   idx = (1:n)';
   idx = idx(:,ones(2.^n,1));
   pranks = sum(allposs.*idx,1);
   tail = 2*length(find(pranks <= w)); % two side.

   % Avoid p>1 if w is in the middle and is double-counted
   p = min(1, tail./(2.^n));
end

if nargout > 1,
   h = (p<=alpha);
   if (nargout > 2)
      stats.signedrank = w;
   end
end

