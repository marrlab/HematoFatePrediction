function yi = interpolateTrajectory(y,method)

yi = y;
yvalid = ~isnan(y);
ygood = find(yvalid);
ybad = find(~yvalid);

y_interp = interp1(ygood,y(yvalid),ybad,method,'extrap');
yi(ybad) = y_interp;

end