function x = constrain(x, loBnd, upBnd)
if x < loBnd, x = loBnd; end
if x > upBnd, x = upBnd; end
