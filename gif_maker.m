% create gif from a video
id = '10';
path = sprintf('newdata/%s/',id);
f = dir(path);
chdir(path)
%x = imresize(imread(f(4).name),0.6);
x = imread(f(4).name);
S = zeros(size(x,1), size(x,2), size(x,3), length(f)-3);
S = uint8(S);
for i=1:length(f)-3
    S(:,:,:,i) = imread(f(i+3).name);
end

filename = sprintf('../%s.gif', id);
for f = 1:size(S,4)
  [SIf,cm] = rgb2ind(S(:,:,:,f),256);
  if f == 1
    imwrite(SIf,cm,filename,'Loop',Inf,'Delay',0);
  elseif f == size(S,4)
    imwrite(SIf,cm, filename,'WriteMode','append','Delay',1.3);
  else
    imwrite(SIf,cm, filename,'WriteMode','append','Delay',0);
  end
end