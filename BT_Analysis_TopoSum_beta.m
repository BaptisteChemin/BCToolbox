size(data)



F = .5:.5:3;
Fb = round((F/header.xstep)+1);

chunk=zeros(size(data,2),1);
for ch =1:size(data,2)
    line = squeeze(data(1,ch,1,1,1,:));
    chunk(ch) = sum(line(Fb));
end

data(1,:,1,1,1,end)=chunk;
plot(squeeze(data(1,:,1,1,1,end)))