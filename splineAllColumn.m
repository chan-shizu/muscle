function yy = splineAllColumn(x,y,xx)

if (size(x,2) == 1 && size(xx,2) == 1)

    [nump,numq] = size(y);
    yy = [];
    for q=1:numq
        tmp = spline(x,y(:,q),xx);
        yy = [yy, tmp];
    end

else
    error('size(x,2) and size(xx,2) must be 1')
end