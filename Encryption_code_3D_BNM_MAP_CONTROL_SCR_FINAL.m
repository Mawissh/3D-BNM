%% Encryption coding:::

imgOrg= double(imread('D:\vivek new article on quantum\3-D chaotic map\Images\Barbara.tif'));

[row,column] = size(imgOrg);
%imgOrg = Img(:,:,1);          % Red   Component
%imgOrg = Img(:,:,2);          % Green Componenet  %%% Color component
%imgOrg = Img(:,:,3);          % Blue Component
%figure(1); imshow(uint8(Img)); title('plaintext image');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chaotic map portion::::::::::
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows, columns, numberOfColorChannels] = size(imgOrg);
maxAllowableSize = 1024;
iteration = 1;
%%Initialize image.
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
    if(iteration==348) 
        break;
    end
end
Arnold_scr= currentScrambledImage;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform the quantum XOR simulation
Diffuse_img = quantumXORSimulation_enc(Img2, Arnold_scr);
con_scr= controlled_qubit_level_scrambling(P, Diffuse_img);
Enc_new = uint8(reshape(con_scr,[512,512]));
Enc_img = uint8(Enc_new);
imwrite(Enc_img,'D:\vivek new article on quantum\bit_planes\Test\barbara_key2cipher_change.png');
subplot(1,3,1); imshow(uint8(imgOrg)); title('original image');
subplot(1,3,2); imshow(uint8(Enc_img)); title('Encrypted image');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
