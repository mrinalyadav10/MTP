clear all;
close all;
tic

% Define initial parameters
xdis = [1, 20];
ydis = [1, 10];
rot = [-20, 0];
n = 10;
numit = 30;
nummut = 1;
f = @image_fit; % Ensure this function is correctly implemented
x = xdis(2) - xdis(1);
y = ydis(2) - ydis(1);
t = rot(2) - rot(1);

% Initialize population
ipopx = rand(1,n) * x + xdis(1);
ipopy = rand(1,n) * y + ydis(1);
ipopt = rand(1,n) * t + rot(1);

% Genetic algorithm loop
for it = 1:numit
    for i = 1:n
        fpop(i) = -f(ipopx(i), ipopy(i), ipopt(i));
    end
    maxf(it) = max(fpop);
    meanf(it) = mean(fpop);
    
    % Fitness normalization
    m = min(fpop);
    fpop = fpop - m;
    cpop(1) = fpop(1);
    for i = 2:n
        cpop(i) = cpop(i-1) + fpop(i);
    end
    
    total_fitness = cpop(n);
    
    % Selection
    for i = 1:n
        p = rand * total_fitness;
        j = find(cpop - p > 0, 1);
        if isempty(j)
            j = n;
        end
        parentx(i) = ipopx(j);
        parenty(i) = ipopy(j);
        parentt(i) = ipopt(j);
    end
    
    % Crossover
    for i = 1:2:n-1
        r = rand;
        ipopx(i) = r * parentx(i) + (1-r) * parentx(i+1);
        ipopx(i+1) = (1-r) * parentx(i) + r * parentx(i+1);
        ipopy(i) = r * parenty(i) + (1-r) * parenty(i+1);
        ipopy(i+1) = (1-r) * parenty(i) + r * parenty(i+1);
        ipopt(i) = r * parentt(i) + (1-r) * parentt(i+1);
        ipopt(i+1) = (1-r) * parentt(i) + r * parentt(i+1);
    end
    
    % Mutation
    for i = 1:nummut
        z = rand;
        ipopx(ceil(z*n)) = xdis(1) + z * x;
        ipopy(ceil(z*n)) = ydis(1) + z * y;
        ipopt(ceil(z*n)) = rot(1) + z * t;
    end
end

% Final fitness evaluation
for i = 1:n
    fpop(i) = -f(ipopx(i), ipopy(i), ipopt(i));
end
ffinal = max(fpop);
j = find(fpop == ffinal);
if isempty(j)
    j = n;
else
    j = j(1);
end

xfinal = ipopx(j);
yfinal = ipopy(j);
tfinal = ipopt(j);

% Load and process images
IM1 = imread('/Users/Mrinal/Desktop/Image.png');
IM2 = imread('/Users/Mrinal/Desktop/Imager.png');
IM1 = double(IM1);
IM2 = double(IM2);
IM3 = double(IM2);  % For display purpose

% Rotate and resize IM2 based on genetic algorithm results
J = imrotate(IM2, tfinal(1), 'bilinear');

% Resizing and matching dimensions for subtraction
outputSize = [size(J, 1), size(J, 2)];
IM1_resized = imresize(IM1, 'OutputSize', outputSize, 'Method', 'bilinear');

% Ensure both images are of the same class for subtraction
IM2_resized = uint8(J);
IM1_resized_uint8 = uint8(IM1_resized);

% Perform subtraction
c3 = imsubtract(IM2_resized, IM1_resized_uint8);

% Error calculation
err = abs(sum(c3(:)) / numel(c3));

% Error calculation
err = abs(sum(c3(:)) / numel(c3));

% Displaying the images
figure;
subplot(1,3,1), imshow(uint8(IM1)), title('Reference Image');
subplot(1,3,2), imshow(uint8(IM3)), title('Input Image');
subplot(1,3,3), imshow(uint8(IM2_resized)), title('Registered and Processed Image');

% Display error
disp(['Error: ', num2str(err)]);

% Measure execution time
executionTime = toc;
disp(['Execution Time: ', num2str(executionTime), ' seconds']);



function f = image_fit(x, y, t)
    IM1 = imread('/Users/Mrinal/Desktop/Imager.png');
    IM2 = imread('/Users/Mrinal/Desktop/Image.png');
    IM1 = double(IM1);
    IM2 = double(IM2);
    
    J = imrotate(double(IM2), t, 'bilinear');
    J = abs(J) .* 255 / max(J(:));
    [n1, n2] = size(IM1);
    [n3, n4] = size(J);
    
    if ((n3 < n1) && (n4 < n2))
        x = 1:n3;
        y = 1:n4;
        IM1 = IM1(x, y);
    elseif ((n3 < n1) && (n4 >= n2))
        x = 1:n3;
        y = 1:n2;
        IM1 = IM1(x, y);
    elseif ((n3 >= n1) && (n4 < n2))
        x = 1:n1;
        y = 1:n4;
        IM1 = IM1(x, y);
    end
    
    [n1, n2] = size(IM1);
    [n3, n4] = size(J);
    
    if x > n3 - n1
        x = n3 - n1 - 1;
        IM1(1:n1, 1:n2) = 255;
    end
    if y > n4 - n2
        y = n4 - n2 - 1;
        IM1(1:n1, 1:n2) = 255;
    end
    if x < 0
        x = 0;
        IM1(1:n1, 1:n2) = 255;
    end
    if y < 0
        y = 0;
        IM1(1:n1, 1:n2) = 255;
    end
    
    if ((x <= n3 - n1) && (y <= n4 - n2))
        xt = 1:n1;
        yt = 1:n2;
        xx = round(xt + x);
        yy = round(yt + y);
        IM2 = round(J(xx, yy));
    elseif ((x <= n3 - n1) && (y > n4 - n2))
        xt = 1:n1;
        yt = 1:n2;
        xx = round(xt + x);
        yy = round(yt);
        IM2 = round(J(xx, yy));
    elseif ((x > n3 - n1) && (y <= n4 - n2))
        xt = 1:n1;
        yt = 1:n2;
        xx = round(xt);
        yy = round(yt + y);
        IM2 = round(J(xx, yy));
    elseif ((x > n3 - n1) && (y > n4 - n2))
        xt = 1:n1;
        yt = 1:n2;
        xx = round(xt);
        yy = round(yt);
        IM2 = round(J(xx, yy));
    end
    
    rows = size(IM1, 1);
    cols = size(IM2, 2);
    N = 256;
    h = zeros(N, N);
    for ii = 1:rows
        for jj = 1:cols
            h(IM1(ii, jj) + 1, IM2(ii, jj) + 1) = h(IM1(ii, jj) + 1, IM2(ii, jj) + 1) + 1;
        end
    end
    
    [r, c] = size(h);
    b = h ./ (r * c);
    y_marg = sum(b);
    x_marg = sum(b');
    Hy = 0;
    for i = 1:c
        if (y_marg(i) == 0)
            % do nothing
        else
            Hy = Hy + -(y_marg(i) * (log2(y_marg(i))));
        end
    end
    Hx = 0;
    for i = 1:r
        if (x_marg(i) == 0)
            % do nothing
        else
            Hx = Hx + -(x_marg(i) * (log2(x_marg(i))));
        end
    end
    h_xy = -sum(sum(b .* (log2(b + (b == 0)))));
    f = -(Hx + Hy - h_xy);
end
