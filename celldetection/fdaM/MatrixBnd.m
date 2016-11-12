function xmat = MatrixBnd(xmat, loBnd, upBnd)
ncol = size(xmat,2);
errwrd = 0;
if length(loBnd) == 1
    loBnd = loBnd*ones(1,ncol);
else
    if length(loBnd) ~= ncol
        errwrd = 1;
    end
end
if length(upBnd) == 1
    upBnd = upBnd*ones(1,ncol);
else
    if length(upBnd) ~= ncol
        errwrd = 1;
    end
end
if errwrd
    error('loBnd or upBnd of wrong length');
end
for j=1:ncol
    if xmat(:,j) < loBnd(j), xmat(:,j) = loBnd(j); end
    if xmat(:,j) > upBnd(j), xmat(:,j) = upBnd(j); end
end
