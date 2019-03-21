val = zeros(256,256,256,3);
for i=1:numel(val); val(i)=i; end
tic; for i=1:1000; slice = val(:,2,:,:); end; toc