function [phdl, thdl] = ...
    phaseplanePlot(fdobj, Lfdobj1, Lfdobj2, plotlab, plotarg)
%  PHASEPLANEPLOT plots a phase/plane diagram
%  Arguments:
%  FDOBJ   ... A univariate functional data object
%  LFDOBJ1 ... Either a nonnegative integer or a linear differential
%              operator argument for variation of the horizonal axis
%  LFDOBJ2 ... Either a nonnegative integer or a linear differential
%              operator argument for variation of the vertical axis
%  PLOTLAB ... Either a cell array of labels for plotted points, or
%              a string array containing labels
%  PLOTARG ... A vector with same length as the number of rows or cells
%              of PLOTLAB indicating independent variable values to be
%              plotted.  By default equally spaced within the range.
%  Returns:  handle objects:
%  PHDL ... handle for the continuous line that is plotted
%  THDL ... handle for the plotting labels
%  LHDL ... handle for the coordinate axes

%  Last modified 30 September 2009

%  set default derivative indices

if nargin < 5,  plotarg = [];          end
if nargin < 4,  plotlab = [];          end
if nargin < 3,  Lfdobj2 = int2Lfd(2);  end
if nargin < 2,  Lfdobj1 = int2Lfd(1);  end
   
basisobj = getbasis(fdobj);

rng     = getbasisrange(basisobj);
nbasis  = getnbasis(basisobj);
nfine   = min(201,10*nbasis+1);
evalarg = linspace(rng(1), rng(2), nfine)';

% Compute points to plot curve

D1fine = eval_fd(evalarg, fdobj, Lfdobj1);
D2fine = eval_fd(evalarg, fdobj, Lfdobj2);

D1rng = [min(D1fine),max(D1fine)];
D2rng = [min(D2fine),max(D2fine)];

if ~isempty(plotlab)
    if iscell(plotlab), plotlab = char(plotlab);  end
    nlab = size(plotlab,1);
    if isempty(plotarg)
        delta   = (rng(2) - rng(1))/(2*nlab);
        plotarg = linspace(rng(1)+delta, rng(2)-delta, nlab)';
    else
        if length(plotarg) ~= nlab
            error('Number of labels and plotting points do not match.')
        end
    end
else
    nlab = 0;
end

% Set up the plot

if nlab == 0
    phdl = plot(D1fine, D2fine, 'b-');
else
    D1vec = eval_fd(plotarg, fdobj, Lfdobj1);
    D2vec = eval_fd(plotarg, fdobj, Lfdobj2);
    phdl  = plot(D1fine, D2fine, 'b-', D1vec, D2vec, 'bo');
    thdl  = text(D1vec, D2vec, plotlab);    
end

if prod(D1rng) < 0 && prod(D2rng) < 0
    lhdl = line([D1rng', zeros(2,1)], [zeros(2,1), D2rng']);
    set(lhdl, 'LineStyle', '--', 'color', 'b')
end

if prod(D1rng) < 0 && prod(D2rng) < 0
    lhdl = line(D1rng,zeros(1,2));
    set(lhdl, 'LineStyle', '--', 'color', 'b')
end

if prod(D1rng) >= 0 && prod(D2rng) < 0
    lhdl = line(zeros(1,2), D2rng);
    set(lhdl, 'LineStyle', '--', 'color', 'b')
end







  

