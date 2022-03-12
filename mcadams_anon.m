% McAdams Anonymization in MATLAB

function [sig_rec, fs] = mcadams_anon(inpath, mcadams)
% mcadams = the McAdams coefficient

winLenMs = 20;
shiftLenMs = 10;
lp_order = 20;

[insig, fs] = audioread(inpath);
insig = insig(:,1);
[len, ch] = size(insig);

winlen = floor(winLenMs * fs * 0.001);
shiftlen = floor(shiftLenMs * fs * 0.001);

% FFT processing parameters
NFFT = 2^(ceil(log2(winlen)));

% analysis and synthesis window which satisfies the constraint
wPR = hann(winlen);
K = sum(wPR)/shiftlen;
win = sqrt(wPR/K);
Nframes = 1+floor((len-winlen)/shiftlen);

% OLA FFT processing
sig_rec = zeros(1,len);

for m = 1:Nframes
	% Indices of the m-th frame
	index = (m-1)*shiftlen+1:min((m-1)*shiftlen+winlen, len);

	% Windowed m-th frame (other than rectangular window)
	frame = insig(index).*win;

	% Get lpc coefficients
	a_lpc = lpc(frame+eps, lp_order);

	% get poles
	[z,p,k] = tf2zpk(1, a_lpc);
	poles = p;

	% index of imaginary poles
	ind_imag = find(imag(poles) ~= 0);

	% index of first imaginary poles
	ind_imag_con = ind_imag(1:2:length(ind_imag));

	% here we define the new angles of the poles, shifted accordingly to the mcadams coefficient
	% values >1 expand the spectrum, while values <1 constract it for angles>1
	% values >1 constract the spectrum, while values <1 expand it for angles<1
	% the choice of this value is strongly linked to the number of lpc coefficients
	% a bigger lpc coefficients number constraints the effect of the coefficient to very small variations
	% a smaller lpc coefficients number allows for a bigger flexibility
	new_angles = angle(poles(ind_imag_con)).^mcadams;

	% make sure new angles stay between 0 and pi
	new_angles(find(new_angles>=pi)) = pi;
	new_angles(find(new_angles<=0)) = 0;

	% copy of the original poles to be adjusted with the new angles
	new_poles = poles;
	for k = 1:length(ind_imag_con)
		% compute new poles with the same magnitued and new angles
		new_poles(ind_imag_con(k)) = abs(poles(ind_imag_con(k)))*exp(j*new_angles(k));

		% applied also to the conjugate pole
        new_poles(ind_imag_con(k)+1) = abs(poles(ind_imag_con(k)+1))*exp(-j*new_angles(k));

	end

	% recover new, modified lpc coefficients
	a_lpc_new = real(poly(new_poles));

	% get residual excitation for reconstruction
	res = filter(a_lpc,1,frame);

	% reconstruct frames with new lpc coefficient
	frame_rec = filter(1, a_lpc_new, res);
	frame_rec = frame_rec.*win;

	outindex = (((m-1)*shiftlen+1):((m-1)*shiftlen+length(frame_rec)));

	% overlap add
	sig_rec(outindex) = sig_rec(outindex) + frame_rec';

end

sig_rec = sig_rec/max(abs(sig_rec));