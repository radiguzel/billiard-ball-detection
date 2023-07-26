clear;
v = VideoReader('data/v2.mp4');

frameSize = round(v.Duration * v.FrameRate);
width = v.Width;
height = v.Height;

video = uint8(zeros(height, width, 3, frameSize));
for i = 1:frameSize
    frame = readFrame(v);
    video(:,:,:,i) = frame;
end
% save a frame as an example
imwrite(video(:,:,:,1),'first_frame.png');