%% Encryption and Decryption coding:::

imgOrg= double(imread('D:\vivek new article on quantum\3-D chaotic map\Images\Barbara.tif'));

[row,column] = size(imgOrg);
%imgOrg = Img(:,:,1);          % Red   Component
%imgOrg = Img(:,:,2);          % Green Componenet  %%% Color component
%imgOrg = Img(:,:,3);          % Blue Component
%figure(1); imshow(uint8(Img)); title('plaintext image');
n = 524288;
cx = zeros(1,n); cy = zeros(1,n); cz = zeros(1,n); 
key=hash(imgOrg,'SHA256');
a1=sum(sum(imgOrg))/(256*row*column);
k4=double(hex2dec(key(1:16))/2^16);
k3=double(hex2dec(key(17:32))/2^16);
k2=double(hex2dec(key(33:48))/2^16);
k1=double(hex2dec(key(49:64))/2^16);
%save('D:\vivek new article on quantum\3-D chaotic map\Extra\Enc_outDNA\7.1.03_info_initial.mat','a1','key','k4','k2','k3','k1');
cx(1) = mod((k1+k4+a1),1);
cy(1)= mod((k2+k4+a1),1);
cz(1) = mod((k3+k4+a1),1);
a = 2.73; b=0.17+10^-15;                    %%a = 3.72; b = 0.17   ,a=2.73             ;cx(1) = 0.21; cy(1) = 0.131; cz(1) = 0.128; 
for j = 1 : n-1
    cx(j+1) = cos(a*cx(j))*sin(1/(cy(j)*(1-cz(j))^2));
    cy(j+1) = sin(a*cx(j)*cy(j)+b*cz(j));
    cz(j+1) = cx(j);
end
% Initialize P sequence
P = [];

% Assuming all sequences are of equal length
seq_length = min([length(cx), length(cy), length(cz)]); % Get the smallest length among x, y, z

k=3;

% Choose every fifth element from each sequence to form P
for i = k:k:seq_length
    P = [P, cx(i), cy(i), cz(i)];
end

P_seq = P(1:262144);
P1_seq = abs(P);
L = mod(ceil(10^15 .* P_seq), 2^8);
Img2 = reshape(L,[512,512]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % tic;
[rows, columns, numberOfColorChannels] = size(imgOrg);
maxAllowableSize = 1024;
iteration = 1;
% Initialize image.
oldScrambledImage = imgOrg;
% The number of iterations needed to restore the image can be shown never to exceed 3N.
Nz = rows;
a=11;
d=19;
while iteration <= 3 * Nz
	% Scramble the image based on the old image.
	for row = 1 : rows % y
		for cols = 1 : columns % x
			c = mod(cols + (a*row), Nz) + 1; % x coordinate
			r = mod((d*cols) +((a*d+1)*row), Nz) + 1; % y coordinate
			% Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
			currentScrambledImage(row, cols, :) = oldScrambledImage(r, c, :);
		end
	end
	% Make the current image the prior/old one so we'll operate on that the next iteration.
	oldScrambledImage = currentScrambledImage;
	% Update the iteration counter.
	iteration = iteration+1;
    if(iteration==59) 
        break;
    end
end
Arnold_scr= currentScrambledImage;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Arnold_scr=enscramble_arnold(imgOrg);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Gray code pre-diffusion on image which  chaotic sequence generated::::
% %scr_img2 = Gray_code(Img2);
% %%%%%%%%%%%%scr_img= Gray_code(imgOrg);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% DNA encoding process on the scrambled image :::
% %Dna_img = DNA_Convert(Arnold_scr,row,column,0.67,0.45,0.51,3.72,0.15,'Encryption');
% % %%%%Dna_img = convert_DNA(Arnold_scr,row,column, 1.2,0.621,0.57,'Encryption');  
% % con_scr= controlled_qubit_level_scrambling(P, Dna_img);
% % cont_scr = uint8(reshape(con_scr,[512,512]));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Perform the quantum XOR simulation
Diffuse_img = quantumXORSimulation_enc(Img2, Arnold_scr);
con_scr= controlled_qubit_level_scrambling(P, Diffuse_img);
Enc_new = uint8(reshape(con_scr,[512,512]));
Enc_img = uint8(Enc_new);
imwrite(Enc_img,'D:\vivek new article on quantum\bit_planes\Test\barbara_key2cipher_change.png');
subplot(1,3,1); imshow(uint8(imgOrg)); title('original image');
subplot(1,3,2); imshow(uint8(Enc_img)); title('Encrypted image');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Decryption code:::

% imgOrg_d = double(imread('D:\vivek new article on quantum\bit_planes\Test\noise\barbara_SP0.4_enc.png'));
% imshow(imgOrg)
% [rows, columns, numberOfColorChannels] = size(imgOrg_d);
% n = 524288;
% cx = zeros(1,n); cy = zeros(1,n); cz = zeros(1,n);
% info_d1 = load("D:/vivek new article on quantum/3-D chaotic map/Extra/Enc_outDNA/barbara_info_initial.mat");
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% For Key sensitivity Analysis::
% % % % % key = 'a5e3b0436c0c059e98be1b5abd724d8f295a3270738370ea6aa8e122c5720411';
% % % % % a1_d=sum(sum(imgOrg))/(256*row*column);
% % % % % k4_d=double(hex2dec(key(1:16))/2^16);
% % % % % k3_d=double(hex2dec(key(17:32))/2^16);
% % % % % k2_d=double(hex2dec(key(33:48))/2^16);
% % % % % k1_d=double(hex2dec(key(49:64))/2^16);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% a1_d=info_d1.a1; k1_d=info_d1.k1;k2_d=info_d1.k2;k3_d=info_d1.k3; k4_d=info_d1.k4;
% cx(1) = mod((k1_d+k4_d+a1_d),1);
% cy(1)= mod((k2_d+k4_d+a1_d),1);
% cz(1) = mod((k3_d+k4_d+a1_d),1);
% 
% a = 2.73; b=0.17;                                 %%a = 3.72; b = 0.19;cx(1) = 0.21; cy(1) = 0.131; cz(1) = 0.128;  
% for j = 1 : n-1
%     cx(j+1) = cos(a*cx(j))*sin(1/(cy(j)*(1-cz(j))^2));
%     cy(j+1) = sin(a*cx(j)*cy(j)+b*cz(j));
%     cz(j+1) = cx(j);
% end
% % Initialize P sequence
% P = [];
% 
% % Assuming all sequences are of equal length
% seq_length = min([length(cx), length(cy), length(cz)]); % Get the smallest length among x, y, z
% 
% k=3;
% 
% % Choose every fifth element from each sequence to form P
% for i = k:k:seq_length
%     P = [P, cx(i), cy(i), cz(i)];
% end
% 
% P_seq = P(1:262144); 
% L = mod(ceil(10^15 .* P_seq), 2^8);
% Img2_d = reshape(L,[512,512]);
% rev_cont_scr = reverse_qubit_level_scrambling(P,imgOrg_d);
% rev_XOR = quantumXORSimulation_dec(rev_cont_scr, Img2_d);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%% 
% %%%%scrd_img2 = Gray_code(Img2_d);
% %%%%rev_Cont_img = reshape(rev_cont_scr,[1,262144]);
% %%%%Dec2_img = arrayfun(@(x) reverseXORshift(x), Dec1_img);
% %%%%%%%%%%Dec3_img =arrayfun(@(x) reverse_qubit_shift_operation(x),Dec2_img);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rev_DNA_img = reverse_qubit_level_scrambling(P,rev_cont_scr);
% rev_Dna_img = reshape(rev_DNA_img,[1,262144]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rev_Dna_convert = DNA_Convert(rev_Dna_img,rows,columns,0.67,0.45,0.51,3.72,0.15,'Decryption');
% %%%%%%%%%%rev_Dna_convert = convert_DNA(rev_Dna_img,row,column,1.2, 0.621,0.57,'Decryption');
% 
%  rev_Arnold_scr = reshape(rev_XOR,[512,512]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%De_arnold_scr=descramble_arnold(rev_Dna);
% 
% maxAllowableSize = 1024;
% iteration = 1;
% % Initialize image.
% oldScrambledImage = rev_Arnold_scr;
% % The number of iterations needed to restore the image can be shown never to exceed 3N.
% Nz = rows;
% a=11;
% d=19;
% while iteration <= 3 * Nz
% 	% Scramble the image based on the old image.
% 	for row = 1 : rows % y
% 		for cols = 1 : columns % x
%             c = mod(cols + (a*row), Nz) + 1; % x coordinate
% 			r = mod((d*cols) +((a*d+1)*row), Nz) + 1; % y coordinate
% 			% % Move the pixel.  Note indexes are (row, column) = (y, x) NOT (x, y)!
% 			currentrevScrambledImage(r, c, :) = oldScrambledImage(row, cols, :);
% 		end
%   end
% 	oldScrambledImage = currentrevScrambledImage;
% 	% Update the iteration counter.
% 	iteration = iteration+1;
%     if(iteration==59)
%         break;
%     end
% end
% De_Arnold_rev_scr= currentrevScrambledImage;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% Decoded_img = uint8(De_Arnold_rev_scr);
% imwrite(Decoded_img,'D:\vivek new article on quantum\bit_planes\Test\noise\dec_barbara_SP0.4.png');
% subplot(1,3,2); imshow(uint8(imgOrg_d)); title('Encrypted image');
% subplot(1,3,3); imshow(uint8(Decoded_img)); title('Reconstruction image');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B_new = quantumXORSimulation_enc(A, B)
    % Ensure A and B are 4x4 and have values within 0-255
    assert(all(size(A) == [512, 512]) && all(size(B) == [512, 512]));
    assert(all(A(:) >= 0 & A(:) <= 255) && all(B(:) >= 0 & B(:) <= 255));

    % Initialize B_new
    B_new = B;

    % Loop through each element in A and B
    for i = 1:512
        for j = 1:512
            % Convert to 8-bit binary
            A_binary = de2bi(A(i, j), 8, 'left-msb');
            B_binary = de2bi(B(i, j), 8, 'left-msb');

            % Perform bitwise quantum XOR simulation
            for k = 1:8
                if A_binary(k) == 1
                    B_binary(k) = ~B_binary(k);
                end
            end

            % Convert B_binary back to decimal and store in B_new
            B_new(i, j) = bi2de(B_binary, 'left-msb');
        end
    end
end

function B_new = quantumXORSimulation_dec(A, B)
    % Ensure A and B are 4x4 and have values within 0-255
    assert(all(size(A) == [512, 512]) && all(size(B) == [512, 512]));
    assert(all(A(:) >= 0 & A(:) <= 255) && all(B(:) >= 0 & B(:) <= 255));

    % Initialize B_new
    B_new = B;

    % Loop through each element in A and B
    for i = 1:512
        for j = 1:512
            % Convert to 8-bit binary
            A_binary = de2bi(A(i, j), 8, 'right-msb');
            B_binary = de2bi(B(i, j), 8, 'right-msb');

            % Perform bitwise quantum XOR simulation
            for k = 1:8
                if A_binary(k) == 1
                    B_binary(k) = ~B_binary(k);
                end
            end

            % Convert B_binary back to decimal and store in B_new
            B_new(i, j) = bi2de(B_binary, 'right-msb');
        end
    end
end