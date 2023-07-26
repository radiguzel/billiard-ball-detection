% This file processes a fixed-camera billiard video and detect the circles. Moreover, it tracks the whiteball and shows its projectory.
% This project is done for the digital image processing course at Bogazici University

clear
% load the video file
load('video.mat')

% table detection PART
mean_im = uint8(mean(video,4));
%mean_im = imresize(mean_im, 0.8);

HSV = rgb2hsv(mean_im);

% find background color in the image
[m,n,~] = size(HSV);
midx = floor(m/2);
midy = floor(n/2);
lenx = floor(m/20);
leny = floor(n/20);
HSV_H= HSV(midx-lenx:midx+lenx, midy-leny:midy+leny,1);
HSV_H_part = reshape(HSV_H, [size(HSV_H,1)*size(HSV_H,2) 1]);
%Find the median Hue
dominant = double(median(HSV_H_part));
hue_threshold = 0.05;
%find the mask
mask = false(m,n);
mask(abs(HSV(:,:,1) - dominant) > hue_threshold) = true;


% find the table corner locations
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

% parameters for detecting the table corner
kup = 0.5;
kdown = 0.8;
kverd = 0.1;
kveru = 0.9;
up = 1;
down = m;
leftd = 1;
leftu = 1;
rightd = n;
rightu = n;
thresholdup = max_hor * kup;
thresholddown = max_hor * kdown;
thresholdverd = max_ver * kverd;
thresholdveru = max_ver * kveru;

% find upperside of the table
for i = 1:m-bin
    total = sum(sum(1 - mask(i:i+bin,:)));
    if total > thresholdup
        up = i;
        break;
    else
        mask(i,:) = 1;
    end
end

% find bottomside of the table 
for i = m:-1:bin+1
    total = sum(sum(1 - mask(i-bin:i,:)));
    if total > thresholddown
        down = i;
        break;
    else
        mask(i,:) = 1;
    end
end

% find leftmost location of the table
for i = 1:n-bin
    total = sum(sum(1 - mask(:,i:i+bin)));
    if total > thresholdverd
        if leftd == 1
            leftd = i;
        end
        if total > thresholdveru
            leftu = i;
            break
        end
    end
end

%find rightmost location of the table
for i = n:-1:bin+1
    total = sum(sum(1 - mask(:,i-bin:i)));
    if total > thresholdverd
        if rightd == size(mask,2)
            rightd = i;
        end
        if total > thresholdveru
            rightu = i;
            break
        end
    end
end

% detect the circles
tic
delay = 0.000001;

% w_centers for keeping the location of white balls
w_centers = [];
total_wrong_ball_detected = 0;

number_of_frames = 0;
hue_threshold = 0.03;
saturation_thre = 0.33;
value_thre = 0.33;
%detecting balls starts here
k_prev = 1;
similar_frames = 0;
frame_increase = 3;
frame_current = frame_increase*k_prev + 1;
frame_fast_count = 0;
white_fast = 0;
white_fast_found = 0;

r = 150;
r_l = 150;
   
% process each frame
while frame_current < size(video, 4)
    colored = video(:,:,:,frame_current);
    hsv = rgb2hsv(colored);
    hsv_h = hsv(:,:,1);
    hsv_h = hsv_h(up:down,leftd:rightd);    
    hsv_s = hsv(:,:,2);
    hsv_s = hsv_s(up:down,leftd:rightd);
    hsv_v = hsv(:,:,3);
    hsv_v = hsv_v(up:down,leftd:rightd);
    
    colored_pr = video(:,:,:,frame_current - frame_increase*k_prev);
    hsv_pr = rgb2hsv(colored_pr);
    hsv_pr_h = hsv_pr(:,:,1);
    hsv_pr_h = hsv_pr_h(up:down,leftd:rightd);
    hsv_pr_s = hsv_pr(:,:,2);
    hsv_pr_s = hsv_pr_s(up:down,leftd:rightd);
    hsv_pr_v = hsv_pr(:,:,3);
    hsv_pr_v = hsv_pr_v(up:down,leftd:rightd);
     
    mask_prev = true(down-up+1,rightd-leftd+1);
    mask_prev(abs(hsv_pr_h(:,:) - dominant) < hue_threshold & hsv_pr_s(:,:) > saturation_thre & hsv_pr_v(:,:) > value_thre ) = false;

    mask = true(down-up+1,rightd-leftd+1);
    mask(abs(hsv_h(:,:) - dominant) < hue_threshold & hsv_v(:,:) > 0.3 & hsv_s(:,:) > 0.33) = false;
    thre_corr = 0.996;
    correlation = corr2(mask, mask_prev);
    if(correlation < thre_corr)
        bw1 = edge(mask,'Sobel');       
        frame_g=rgb2ycbcr(colored);
        bw = edge(frame_g(:,:,1),'Sobel');   

        % edge of the table
        bw1 = bw(up:down,leftd:rightd);
        [centers1, radii1] = imfindcircles(bw1, [9 20]);%tune
        if ~isempty(centers1)
            centers = centers1;
            radii = radii1;
            centers = centers+[leftd up];
        end

        % filter circles around the vertical lines
        slopel = (up-down) / (leftu-leftd);
        sloper = (up-down) / (rightu-rightd);
        centersl = centers - [leftd down];
        centersr = centers - [rightd down];
        slope_l = centersl(:,2)./centersl(:,1);
        slope_r = centersr(:,2)./centersr(:,1);
        centers_f = (slope_l > slopel & slope_r < sloper);
        centers = centers(centers_f==1,:);
        radii = radii(centers_f==1);

        % viscircles(centers, radii,'EdgeColor','b');
        % find the white ball
        wlim = 160;
        thre = 700;
        c_r = round(centers);
        r_r = radii;
        white_index = -1;
        white = [255, 255, 255];
        temp = 10000;
        thre_h = 0.02;
        thre_v = 0.4;
        centers_b=[];
        for ii = 1:size(c_r,1)
            x = c_r(ii,1);
            l = floor(radii(ii)/2);
            y = c_r(ii,2);
            cut = colored(y-l:y,x-l:x+l,:);
            mu = mean(cut,1);
            mu = squeeze(mean(mu,2));
            dif = norm(white - mu);
            if dif < temp 
                temp = dif;
                index = ii;
                mu_f = mu;
            end
            mu_hsv = rgb2hsv( im2double(uint8(round(mu))) );
            mu_h = mu_hsv(1);
            dif_b = norm(dominant-mu_h);
            if dif_b <thre_h %%&& mu_hsv(3) > thre_v
                centers_b=[centers_b ii];
            end
        end
        centers(centers_b,:)=[];
        radii(centers_b)=[];

        % c_r(centers_b,:)=[];
        total_wrong_ball_detected = total_wrong_ball_detected + length(centers_b);
        %index = index - length(centers_b(centers_b<index));
        viscircles(centers, radii,'EdgeColor','r');
        % viscircles(centers(centers_b,:), radii(centers_b),'EdgeColor', 'b');
        % color white ball with 
        white_threshold = thre;
        wlim1 = 110;
        
        if norm(white - mu_f) < thre && mu_f(1) > wlim && mu_f(2) > wlim && mu_f(3) > wlim
            white_index = index;
            w_center = c_r(white_index,:);
            w_centers = [w_centers ; w_center];
            w_radii = r_r(white_index);
            viscircles(w_center, w_radii,'EdgeColor','g');

        elseif exist('w_center','var')
            mask_white = false(m,n);
            mask_white(double(colored(:,:,1)) > wlim1 & double(colored(:,:,2)) > wlim1 & double(colored(:,:,3)) > wlim1) = true;
            mask_white_box = mask_white(max(w_centers(end,2) - r,1) : min(w_center(end,2) + r,m) ,  max(w_center(end,1) - r,1) : min(w_center(end,1) + r,n));
            s = regionprops(mask_white_box,'centroid','MajorAxisLength','MinorAxisLength','Orientation');
              regionprop_plot(s,mask_white_box,3);
            se = strel('disk',5);
            mask_white_box_opened = imopen(mask_white_box,se);
            s1 = regionprops(mask_white_box_opened,'centroid','MajorAxisLength','MinorAxisLength','Orientation');
              regionprop_plot(s1,mask_white_box_opened,4);
            white_ind = -1;
            white_ind = find_white_center(s,s1,w_radii);
            figure(1)
            text( n-160, 160, '??','FontSize', 40,'Color','white');
            if white_ind ~= -1
                x = round(s1(white_ind).Centroid(1));
                y = round(s1(white_ind).Centroid(2));
                w_center = [max(w_centers(end,1) - r,1)+x, max(w_center(end,2)-r,1) + y];
                w_centers = [w_centers ; w_center];
                white_fast_found = white_fast_found + 1;
            end
            white_fast = white_fast + 1;
            frame_fast_count = 0;
        end


        % draw ROI lines
        line([leftd leftu], [down up], 'Marker', '.', 'color','g', 'LineWidth', 1.8);
        line([rightd rightu], [down up], 'Marker', '.', 'color','g', 'LineWidth', 1.8);
        line([rightd leftd], [down down], 'Marker', '.', 'color','g', 'LineWidth', 1.8);
        line([leftu rightu], [up up], 'Marker', '.', 'color','g', 'LineWidth', 1.8);
        if length(w_centers)>2
            line([w_centers(end,1) w_centers(end-1,1)],[w_centers(end,2) w_centers(end-1,2)],'Marker', '.', 'Color','w', 'LineWidth', 4);
        end
        number_of_frames = number_of_frames + 1;
        pause(delay);
        
        frame_fast_count = frame_fast_count + 1;
        if frame_fast_count < 3
            frame_increase = 1;
        else
            frame_increase = 3;
        end
        frame_current = frame_current + frame_increase;
    end
end

white_ball_find_percentage = length(w_centers) / number_of_frames * 100;
white_fast_found_percentage = white_fast_found / white_fast * 100;
similar_frames_percentage = similar_frames / number_of_frames * 100;
dist_thre = 350;

% plot the white line route
% wi_centers1 = smooth(w_centers(:,1));
% wi_centers2 = smooth(w_centers(:,2));
% wi_centers = [wi_centers1 wi_centers2];
wi_centers = w_centers;
for j = 2:length(wi_centers)
     if norm(wi_centers(j,:)-wi_centers(j-1,:)) < dist_thre
        line([wi_centers(j,1) wi_centers(j-1,1)],[wi_centers(j,2) wi_centers(j-1,2)],'Marker', '.', 'Color','w', 'LineWidth', 4);
     end
end
toc 


function regionprop_plot(s, mask, fig)
    figure(fig)
    imshow(mask,'InitialMagnification','fit')
    t = linspace(0, 2*pi, 50);
    hold on
    for k = 1:length(s)
        a = s(k).MajorAxisLength/2;
        b = s(k).MinorAxisLength/2;
        Xc = s(k).Centroid(1);
        Yc = s(k).Centroid(2);
        phi = deg2rad(-s(k).Orientation);
        x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
        y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
        plot(x,y,'r','Linewidth',5)
    end
    hold off
end
