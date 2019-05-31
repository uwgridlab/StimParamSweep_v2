function y = Weibull(p,x)
%y = Weibull(p,x)
%
%Parameters:  p.b slope
%             p.t threshold yeilding ~80% correct
%             p.g guess rate (50% by default)
%             p.e performance at threshold (.79 by default)
%             x   intensity values.

if ~isfield(p,'g')
    p.g = 0.5;  %chance performance
end

if ~isfield(p,'e')
    p.e = (.5)^(1/3);  %threshold performance ( ~80%)
end

%here it is.
k = (-log( (1-p.e)/(1-p.g)))^(1/p.b);
y = 1- (1-p.g)*exp(- (k*x/p.t).^p.b);

