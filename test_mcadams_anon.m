% Test program for McAdams Anonymization

[outsig, fs] = mcadams_anon('orgsound_001.wav',0.9);
audiowrite('anonsig.wav',outsig,fs);