function spm_write_image(h, V, fname)
%  Write spm image

h.fname = fname;
h.private.dat.fname = h.fname;
spm_write_vol(h, V);