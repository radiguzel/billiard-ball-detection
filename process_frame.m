# This file process a single frame of a fixed-camera billiard video and detects the table corners and the balls

colored = imread('first_frame.png');
I = double(colored);
J = rgb2gray(colored);
J = double(J);
mu = mean(J(:));
[m,n] = size(J);
dist = 30;

K = reshape(colored, [m*n  3]);
dominant = double(mode(K,1));
mask = false(m,n);
for i = 1:m
    for j = 1:n
        dist = norm(reshape(I(i,j,:),1,3) - dominant);
        if(dist > 90 )
            mask(i,j)=true;
        end
    end
end
imshow(mask)

% table corner detection
bin = 4;
max_hor = 0;
for i = 1:m-bin
    total = sum(sum(1 - mask(i:i+bin,:)));
    if total > max_hor
        max_hor = total;
    end
end
max_ver = 0;
for i = 1:n-bin
    total = sum(sum(1 - mask(:,i:i+bin)));
    if total > max_ver
        max_ver = total;
    end
end
kup = 0.5;
kdown = 0.8;
kverd = 0.3;
kveru = 0.95;
up = 1;
down = size(mask,1);
left = 1;
leftc = 1;
right = size(mask,2);
rightc = size(mask,2);
thresholdup = round(max_hor*kup);
thresholddown = round(max_hor*kdown);
thresholdverd = round(max_ver*kverd);
thresholdveru = round(max_ver*kveru);


for i = 1:m-bin
    total = sum(sum(1 - mask(i:i+bin,:)));
    if total > thresholdup
        up = i;
        break;
    end
end

for i = m:-1:bin+1
    total = sum(sum(1 - mask(i-bin:i,:)));
    if total > thresholddown
        down = i;
        break;
    end
end

for i = 1:n-bin
    total = sum(sum(1 - mask(:,i:i+bin)));
    if total > thresholdverd
        if left == 1
            left = i;
        end
        if total > thresholdveru
            leftc = i;
            break
        end
    end
end

for i = n:-1:bin+1
    total = sum(sum(1 - mask(:,i-bin:i)));
    if total > thresholdverd
        if right == size(mask,2)
            right = i;
        end
        if total > thresholdveru
            rightc = i;
            break
        end
    end
end


%% apply structure extraction
bw = edge(mask,'Canny');
figure
imshow(colored)
figure
imshow(bw)

% find circles in the table
[centers,radii] = imfindcircles(bw, [5 15] );
center_f = (centers(:,2)>up & centers(:,2)<down & centers(:,1)>left & centers(:,1)<right ); 
centers_filtered= centers(find(center_f==1),:);
radii_filtered = radii(find(center_f==1));
viscircles(centers, radii,'EdgeColor','b');
viscircles(centers_filtered, radii_filtered,'EdgeColor','r');
%line([down left] , [down right], 'color','g') 
