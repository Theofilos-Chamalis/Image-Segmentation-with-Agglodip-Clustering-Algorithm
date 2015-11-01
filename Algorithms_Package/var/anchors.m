%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% choose k centers by the furthest-first method
% source: fastkmeans package from Matlab File Exchange
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [centers, mincenter, mindist, lower, computed] = anchors(firstcenter,k,data)

    [n,dim]  = size(data);
    centers  = zeros(k,dim);
    lower    = zeros(n,k);
    mindist  = Inf*ones(n,1);
    mincenter= ones(n,1);
    computed = 0;
    centdist = zeros(k,k);

    newcenter= firstcenter;

    for j = 1:k
        fprintf('center: %g\n', j);
        if j > 1,
            [~,i] = max(mindist);
            newcenter = data(i,:);
        end

        centers(j,:) = newcenter;
        centdist(1:j-1,j) = calcdist(centers(1:j-1,:),newcenter);
        centdist(j,1:j-1) = centdist(1:j-1,j)';
        computed = computed + j-1;

        inplay = find(mindist > centdist(mincenter,j)/2);

        newdist = sum(data(inplay,:).^2, 2) - 2*data(inplay,:)*newcenter' + newcenter*newcenter';
        computed = computed + size(inplay,1);
        lower(inplay,j) = newdist;

        move = find(newdist < mindist(inplay));
        shift = inplay(move);
        mincenter(shift) = j;
        mindist(shift) = newdist(move);
    end
end