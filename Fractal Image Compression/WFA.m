function wfa_main
    % Define a test matrix F (should be a power-of-2 square matrix)
    F = magic(4);  % Example matrix, change as needed

    % Call the wfa_encoder function with the test matrix
    [W, FF, bits] = wfa_encoder(F);

    % Display the results
    disp('W:');
    disp(W);
    disp('FF:');
    disp(FF);
    disp('Bits:');
    disp(bits);

    % Decode the WFA to reconstruct the matrix
    n = log2(size(F, 1));
    F_reconstructed = wfa_decoder(W, FF, n);

    % Display the reconstructed matrix
    disp('Reconstructed F:');
    disp(F_reconstructed);
end

function leeves = get_leaves(A)
    n = length(A(1,:));
    if n == 2
        leeves = [A(2,1) A(1,1) A(2,2) A(1,2)];
    else
        leeves = [];
        leeves = [leeves get_leaves(A((n/2 + 1):end, 1:n/2))];
        leeves = [leeves get_leaves(A(1:n/2, 1:n/2))];
        leeves = [leeves get_leaves(A((n/2 + 1):end, (n/2 + 1):end))];
        leeves = [leeves get_leaves(A(1:n/2, (n/2 + 1):end))];
    end
end

function f = downsample(x, n)
    len = length(x);
    f = x;
    for k = 1:n
        x = f;
        f = [];
        for m = 1:length(x)/4
            f(m) = mean(x(4*(m-1)+1:4*(m-1)+4));
        end
    end
end

function f = upsample(x, n)
    if n == 0
        f = x';
    elseif n < 0
        f = downsample(x, -n);
    else
        len = length(x);
        for m = 0:len - 1
            f((4^n)*m+1:(4^n)*m+4^n) = x(m+1);
        end
    end
end

function [W, FF, bits] = wfa_encoder(F)
    % Recursive algorithm with cosine initial basis
    disp('Running WFA Encoder...');
    global W FF n depth G tree state SymTable SymCount init_states init_res;
    init_states = 256;
    init_res = 4;
    W = [];
    W(:,:,1) = zeros(1, init_states + 1); 
    W(:,:,2) = zeros(1, init_states + 1);
    W(:,:,3) = zeros(1, init_states + 1); 
    W(:,:,4) = zeros(1, init_states + 1);
    FF = [];
    for kk = 1:init_states
        FF = [FF ; mean(400*cos((kk-1)*linspace(0, pi, 256))')];
    end
    FF = [FF ; mean(mean(F))];
    SymTable = [0 ; 4*(init_states + 1)];
    SymCount = 4*(init_states + 1);
    G = 180;
    depth = log2(length(F(1,:)));
    tree{depth + 1} = get_leaves(F);
    for level = depth:-1:1
        tree{level} = downsample(tree{level + 1}, 1);
    end
    n = 1;
    state = [1; 1];
    totalcost = make_wfa(1, depth, Inf);
    % Check for no zeros
    if SymTable(2, 1) == 0
        bits = SymTable(2, 2:end) * log2(SymCount ./ SymTable(2, 2:end))';
    else
        bits = SymTable(2,:) * log2(SymCount ./ SymTable(2,:))';
    end
    bits = bits + 16 * n;
end

function cost = make_wfa(ii, k, maks)
    global W FF n depth G tree state SymTable SymCount init_states init_res;
    if maks <= 0 || k == 0
        cost = Inf; return
    end
    cost = 0;
    f = [];
    for kk = 1:init_states
        f = [f upsample(400*cos((kk-1)*linspace(0, pi, 256))', k - init_res - 1)'];
    end
    for a = 0:3
        pos = state(:, ii);
        pos_psi(1) = pos(1) + 1;
        pos_psi(2) = 4*(pos(2) - 1) + 1 + a;
        psi = tree{pos_psi(1) + k - 1}(((4^(k-1))*(pos_psi(2)-1) + 1):((4^(k-1))*(pos_psi(2)-1) + 1 + 4^(k-1) - 1))';
        PHI = [];
        for m = 1:n
            pos = state(:, m);
            if pos(1) + k - 1 <= depth + 1
                PHI = [PHI ; tree{pos(1) + k - 1}(((4^(k-1))*(pos(2)-1) + 1):((4^(k-1))*(pos(2)-1) + 1 + 4^(k-1) - 1))];
            else
                overshoot = pos(1) + k - 1 - (depth + 1);
                temp1 = tree{depth + 1}(((4^(k-1-overshoot))*(pos(2)-1) + 1):((4^(k-1-overshoot))*(pos(2)-1) + 1 + 4^(k-1-overshoot) - 1));
                PHI = [PHI ; upsample(temp1, overshoot)];
            end
        end
        PHI = [f PHI'];
        r = PHI \ psi;
        r = r';
        % Quantize
        r = double(int32(r * (2^(k + 4)))) / (2^(k + 4));
        s1 = SymTable(2,:) * log2(SymCount ./ SymTable(2,:))';
        NewSymTable = SymTable;
        NewSymTable(2, 1) = NewSymTable(2, 1) - length(r);
        file = unique(r);
        for p = 1:length(file)
            tablepos = find(NewSymTable(1,:) == file(p));
            occur = length(find(r == file(p)));
            if isempty(tablepos)
                NewSymTable = [NewSymTable [file(p) ; occur]];
            else
                NewSymTable(2, tablepos) = NewSymTable(2, tablepos) + occur;
            end
        end
        % Check for no zeros
        if NewSymTable(2, 1) == 0
            s2 = NewSymTable(2, 2:end) * log2(SymCount ./ NewSymTable(2, 2:end))';
        else
            s2 = NewSymTable(2,:) * log2(SymCount ./ NewSymTable(2,:))';
        end
        s = s2 - s1;
        cost1 = G * s + norm(psi - PHI * r', 2)^2;
        n0 = n;
        n = n + 1;
        state = [state pos_psi'];
        FF = [FF ; tree{pos_psi(1)}(pos_psi(2))];
        W(1, n + init_states, 1) = 0;
        W(n, 1, 1) = 0;
        W(ii, n + init_states, a + 1) = 1;
        OldSymCount = SymCount;
        SymCount = SymCount + 4 * (n0 + init_states) + 4 * (n);
        tablepos = find(SymTable(1,:) == 1);
        if isempty(tablepos)
            SymTable = [SymTable [1 ; 1]];
        else
            SymTable(2, tablepos) = SymTable(2, tablepos) + 1;
        end
        SymTable(2,1) = SymTable(2,1) + 4 * (n0 + init_states) + 4 * (n) - 1;
        s2 = SymTable(2,:) * log2(SymCount ./ SymTable(2,:))';
        s = s2 - s1 + 16;
        cost2 = G * s + make_wfa(n, k - 1, min([cost1 maks - cost]));
        cost = cost + min([cost1 cost2]);
        if cost == cost2
            state = state(:, 1:n0);
            SymTable = NewSymTable;
            SymCount = OldSymCount;
            W = W(1:n0 + init_states, 1:n0 + init_states, :);
        end
        n = n0;
    end
    if cost > maks
        cost = Inf;
    end
end

function F = wfa_decoder(W, FF, n)
    disp('Running WFA Decoder...');
    global F pos leeves;
    init_states = 256;
    init_res = 4;
    total_states = length(FF);
    level0 = FF';
    
    % Initialize level1 using downsampled initial states
    level{1} = [];
    for kk = 1:init_states
        level{1} = [level{1} downsample(400*cos((kk-1)*linspace(0, pi, 256))', init_res - 1)];
    end
    
    % Construct further levels
    for k = 2:(n-1)
        level{k} = zeros(1, (4^k) * total_states);
        assign = [];
        for kk = 1:init_states
            assign = [assign upsample(400*cos((kk-1)*linspace(0, pi, 256))', k - init_res)];
        end
        level{k}(1:init_states*(4^k)) = assign;
    end
    
    % Decode using WFA coefficients
    for k = 2:(n-1)
        for state = (init_states + 1):total_states
            pos = 0;
            aaa = 0;
            for iterate = 1:4^(k-2)
                for a = 0:3
                    aa = 0;
                    for i = ((state-1)*(4^k) + 1 : 4^(k-1) : state*(4^k)) + a + pos
                        thingy = level{k-1}((1 : 4^(k-1) : end) + aaa);
                        if length(W(state - init_states, :, aa + 1)) == length(thingy)
                            leeves(i) = W(state - init_states, :, aa + 1) .* thingy';
                        end
                        aa = mod(aa + 1, 4);
                        if mod(aa, 4) == 0
                            aaa = aaa + 1;
                        end
                    end
                end
                pos = 4 * iterate;
            end
        end
    end
    
    % Final level
    k = n;
    leeves = zeros(1, 4^k);
    state = init_states + 1;
    pos = 0;
    aaa = 0;
    for iterate = 1:4^(k-2)
        for a = 0:3
            aa = 0;
            for i = (1 : 4^(k-1) : 4^k) + a + pos
                thingy = level{k-1}((1 : 4^(k-1) : end) + aaa);
                if length(W(state - init_states, :, aa + 1)) == length(thingy)
                    leeves(i) = W(state - init_states, :, aa + 1) .* thingy';
                end
                aa = mod(aa + 1, 4);
                if mod(aa, 4) == 0
                    aaa = aaa + 1;
                end
            end
        end
        pos = 4 * iterate;
    end
    
    % Pixel reconstruction
    pos = 1;
    F = zeros(2^n, 2^n);
    pixel(n, 1, 1);
    
    % Recursive function to place leaves into F
    function [] = pixel(nn, row, col)
        if nn == 1
            F(row + 1, col) = leeves(pos);
            F(row, col) = leeves(pos + 1);
            F(row + 1, col + 1) = leeves(pos + 2);
            F(row, col + 1) = leeves(pos + 3);
            pos = pos + 4;
        else
            pixel(nn - 1, row + 2^(nn-1), col);
            pixel(nn - 1, row, col);
            pixel(nn - 1, row + 2^(nn-1), col + 2^(nn-1));
            pixel(nn - 1, row, col + 2^(nn-1));
        end
    end
end


