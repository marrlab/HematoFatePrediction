function [yi,SSE,xfd,lambda,outliers] = smoothTrajectoryWoLambda(xi,y_raw)

assert(numel(xi) == numel(y_raw))

bfuns = round(numel(xi)/5);

bspl_basis = create_bspline_basis([xi(1) xi(end)],bfuns,4);

%  estimate standard error of fit
loglam = (4:.1:6)';
nlam   = length(loglam);

dfsave  = zeros(nlam,1);
gcvsave = zeros(nlam,1);

%  loop through smoothing parameters

for ilam=1:nlam
    lambda = 10^loglam(ilam);
%     display(['lambda = ',num2str(lambda)])
    fdParobj = fdPar(bspl_basis, 2, lambda);
    [xfd, df, gcv] = smooth_basis(xi, y_raw, fdParobj);
    dfsave(ilam)  = df;
    gcvsave(ilam) = sum(gcv);
    disp(lambda)
%     plotfit_fd(y_raw, xi, xfd)
end

% doplot = 1;
% if doplot
%     disp('Log lambda    df          gcv')
% disp([loglam, dfsave, gcvsave])
% 
% subplot(2,1,1)
% plot(loglam, gcvsave, 'o-')
% ylabel('\fontsize{16} GCV Criterion')
% title('\fontsize{16} Temperature Smoothing')
% subplot(2,1,2)
% plot(loglam, dfsave, 'o-')
% xlabel('\fontsize{16} log_{10} \lambda')
% ylabel('\fontsize{16} Degrees of Freedom')
% end
[~,minidx] = min(gcvsave);
lambda = 10^loglam(minidx);
par = fdPar(bspl_basis,2,lambda);
[xfd, df, gcv, coef, SSE] = smooth_basis(xi,y_raw,par);
yi_fitted = eval_fd(xi,xfd);
residual_diffs = y_raw-yi_fitted;

ubound = quantile(residual_diffs,.98);
lbound =  quantile(residual_diffs,.08);

outliers = residual_diffs < lbound | residual_diffs > ubound;
yi = getcoef(xfd)';

% plot(residual_diffs,'.')
% hold on
% line([min(xi),max(xi)],[ubound ubound],'Color','r')
% line([min(xi),max(xi)],[lbound lbound],'Color','r')



end