nDims = 3;
nGPoints = [256 256 256];
nGPointsProd = [1 cumprod(nGPoints)];

val = zeros(nGPoints);
val(:) = 1:prod(nGPoints);

nInt = prod(nGPoints-2);
nExt = prod(nGPoints)-nInt;

N = 1000;

tic
for i=0:N

    v = val(1,:,:);
    halo = reshape(v,1,numel(v));
    v = val(end,:,:);
    halo = [halo reshape(v,1,numel(v))];
    v = val(2:end-1,1,:);
    halo = [halo  reshape(v,1,numel(v))];
    v = val(2:end-1,end,:);
    halo = [halo  reshape(v,1,numel(v))];
    v = val(2:end-1,2:end-1,1);
    halo = [halo  reshape(v,1,numel(v))];
    v = val(2:end-1,2:end-1,end);
    halo = [halo  reshape(v,1,numel(v))];
    
end
toc