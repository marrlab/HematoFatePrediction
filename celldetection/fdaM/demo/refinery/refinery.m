%  -------------------------------------------------------------
%            Analyses of the Corpus Christi refinery data
%  -------------------------------------------------------------

%  input the data

load refinery.dat

t = refinery(:,1);  %  observation time
u = refinery(:,2);  %  reflux flow
x = refinery(:,3);  %  tray 47 level

%  plot data

subplot(2,1,1)
plot(t, x, 'k.')
ylabel('\fontsize{12} Tray 47 level')
axis([0,193,215,220])

subplot(2,1,2)
plot(t, u, 'k.')
xlabel('\fontsize{12} Time (min)')
ylabel('\fontsize{12} Reflux flow')
axis([0,193,19.8,20.6])

% print -dps2 'c:\MyFiles\fdabook\revision\figs.dir\refinerydata.ps'

% center the data

u = u - mean(u(1:60));
x = x - mean(x(1:60));

n = length(t);

range = [t(1), t(n)];
tval = t;
delta = 1/(n-1);
trng = [t(1), t(n)];

tbreak = t(67);

%  set up basis for input variable

norder = 1;
nubasis = 2;
ubreaks = [trng(1), tbreak, trng(2)];
ubasis = create_bspline_basis(trng, nubasis, norder, ubreaks);

ufd = data2fd(u, tval, ubasis);

%  set up basis for the output variable
%  put three coincident knots at tbreak

norder  =  4;
xknots  = [linspace(0, tbreak, 3), tbreak, linspace(tbreak, trng(2), 5)];
nxbasis = length(xknots) + norder - 2;
xbasis  = create_bspline_basis(trng, nxbasis, norder, xknots);

xfd = data2fd(x, tval, xbasis);
xvec = eval_fd(tval, xfd);

%  plot the data with fit and knot locations

ahdl = axes('Box', 'on', 'FontSize', 13);
set(ahdl, 'Xlim', [0,193]);
set(ahdl, 'Ylim', [-0.5,4.5]);
set(ahdl, 'Xtick', 0:50:150);
set(ahdl, 'Ytick', 0:4);
lhdl = plot(tval, x, 'b*', tval, xvec, 'r-');
set(lhdl, 'LineWidth', 2)
xlabel('Time (minutes)', 'FontSize', 16);
ylabel('Tray 47 level',  'FontSize', 16);
for i=2:length(xknots)-1
    lhdl = line([xknots(i),xknots(i)], [-0.5,4.5]);
    set(lhdl, 'LineWidth', 1, 'LineStyle', '--', 'color', 'g')
end

print -dps2 'c:/MyFiles/fdabook/revision/figs.dir/refinery.ps'

