%-----------------------------------------------------------------------
% Job saved on 16-Nov-2022
% spm SPM - SPM12 (7487)
% SB 
%-----------------------------------------------------------------------
global Trk
global parentf
global atlas

matlabbatch{1}.spm.util.imcalc.input = {
                                        [atlas,',1']
                                        [parentf,Trk,',1']
                                        };
matlabbatch{1}.spm.util.imcalc.output = ['m_',Trk];
matlabbatch{1}.spm.util.imcalc.outdir = {'E:\Research_Projects\MEG_CamCAN\_Glasser52_to_CABNP\3_Multiplied'};
matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
matlabbatch{1}.spm.util.imcalc.options.dtype = 8;


