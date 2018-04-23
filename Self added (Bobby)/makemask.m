function mask = makemask(FOVx,FOVy)

mask = zeros(FOVx,FOVy);
xcenter = FOVx/2;
ycenter = FOVy/2;
dists = zeros(FOVx,FOVy);
for ii = 1:FOVx;
    for jj = 1:FOVy;
        dists(ii,jj) = sqrt((ii-xcenter)^2.0 + (jj - ycenter)^2.0);
    end
end

maxdist = max(dists(:));

for ii = 1:FOVx;
    for jj = 1:FOVy;
        p = -dists(ii,jj)/(0.5*maxdist) + 1; % For linear probability density
%         p = exp(-(-log(0.2)/maxdist)*dists(ii,jj)); % For exponentially decaying probability density.
        randomnumber = rand;
        if randomnumber < p
            mask(ii,jj) = 1;
        else
            mask(ii,jj) = 0;
        end
    end
end