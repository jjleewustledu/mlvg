%-- 2/4/24, 23:19 --%
fn
isfile(fn)
sum(p.select_vec)
208*300*320
sum(p.select_vec)^(1/3)
440*440*165
ic = mlfourd.ImagingContext2(t1w_fqfn)
ic.view
med
wmparc_on_obj
obj
med
med.bids
med.imagingContext
bk
wmparc_on_obj
out
ic1.filepath = targ_pth
ic1.save
ic1
close
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108293/ses-20210421152358/pet')
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421152358_trc-ho_proc-BrainMoCo2-createNiftiMovingAvgFrames-ParcWmparc-reshape-to-wmparc-select-all.nii.gz")
ic.addJsonMetadata(struct("martinv1", 0.05))
ic.json_metadata
ic.save
ic = mlfourd.ImagingFormatContext2("sub-108293_ses-20210421152358_trc-ho_proc-BrainMoCo2-createNiftiMovingAvgFrames-ParcWmparc-reshape-to-wmparc-select-all.nii.gz")
ic.addJsonMetadata(struct("martinv1", 0.05))
ic.json_metadata
ic.save
ic = mlfourd.ImagingFormatContext2("sub-108293_ses-20210421152358_trc-ho_proc-BrainMoCo2-createNiftiMovingAvgFrames-ParcWmparc-reshape-to-wmparc-wmparc.nii.gz")
ic.addJsonMetadata(struct("martinv1", 0.05))
ic.json_metadata
ic.save
open mlpet.Radionuclides
close all
.011*60
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421152358_trc-ho_proc-BrainMoCo2-createNiftiMovingAvgFrames-ParcWmparc-reshape-to-wmparc-wmparc.nii.gz")
ic = mlfourd.ImagingFormatContext2("sub-108293_ses-20210421152358_trc-ho_proc-BrainMoCo2-createNiftiMovingAvgFrames-ParcWmparc-reshape-to-wmparc-wmparc.nii.gz")
ic.img = ic.img(1:3,:);
ic.fileprefix = "sub-108293_ses-20210421152358_trc-ho_proc-BrainMoCo2-createNiftiMovingAvgFrames-ParcWmparc-reshape-to-wmparc-only3"
ic.save
plot(ic.img)
plot(ic.img')
ifc = mlfourd.ImagingFormatContext2("sub-108293_ses-20210421152358_trc-ho_proc-BrainMoCo2-createNiftiMovingAvgFrames-ParcWmparc-reshape-to-wmparc-only3_dynesty-Raichle1983Model-rho-pred.nii.gz")
figure; plot(ifc.img')
close all
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108293/ses-20210421155709/pet')
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421155709_trc-fdg_proc-MipIdif_idif.nii.gz")
plot(ic)
close
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421155709_trc-fdg_proc-MipIdif_idif_dynesty-Boxcar-ideal.nii.gz")
plot(ic)
ic1 = mlfourd.ImagingContext2("sub-108293_ses-20210421155709_trc-fdg_proc-MipIdif_idif.nii.gz")
plot(ic1)
close all
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108237/ses-20221031100910/pet')
ic = mlfourd.ImagingContext2("sub-108237_ses-20221031100910_trc-co_proc-MipIdif_idif_dynesty-Boxcar-ideal.nii.gz")
plot(ic)
close all
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421144815_trc-co_proc-MipIdif_idif_dynesty-Boxcar-ideal.nii.gz")
plot(ic)
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108293/ses-20210421144815/pet')
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421144815_trc-co_proc-MipIdif_idif_dynesty-Boxcar-ideal.nii.gz")
plot(ic)
close
close all
558*60
558*60*1.05
1.05*558/60
ic
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421144815_trc-co_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic)
t_drawn = []; activity = []
t_drawn = {}; activity = []
t_drawn = []; activity = []
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2021,4,21))
rm.fdg.TIMEDRAWN_Hh_mm_ss
rm.countsFdg.TIMEDRAWN_Hh_mm_ss
t_drawn = rm.countsFdg.TIMEDRAWN_Hh_mm_ss;
activity = rm.countsFdg.DECAYCorrSpecificActivity_KBq_mL;
T = table(t_drawn, activity)
T(7:8,:) = [];
T([2,4,6,end],:) = [];
rm.tracerAdmin
rm.tracerAdmin{5,2}
t_elapsed = seconds(t_drawn - rm.tracerAdmin{5,2})
T = table(t_drawn, activity)
T = table(t_drawn, t_elapsed, activity)
T(7:8,:) = [];
T([2,4,6,end],:) = [];
ic.json_metadata.times(end)
ic.json_metadata.timesEnd(end)
ic.json_metadata.timesMid(end)
ic.json_metadata
ic.json_metadata.timesMid(1)
ic.json_metadata.times(1)
close all
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421155709_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
ic.json_metadata
ic.json_metadata.timesMid(end)
ic.json_metadata.times(end)
ifc = ic.imagingFormat;
ifc.img(end-7:end)
plot(ic)
ifc.img(end-17:end)
ifc.json_metadata.timesMid(1)
t_elapsed = seconds(t_drawn - datetime(20214,21,15,57,9)) + 31
t_elapsed = seconds(t_drawn - datetime(2021,4,21,15,57,9)) + 31
t_elapsed = seconds(t_drawn - datetime(2021,4,21,15,57,9, TimeZone="local")) + 31
T = table(t_drawn, t_elapsed, activity)
T(7:8,:) = [];
T([2,4,6,end],:) = [];
T.activity = 1e3*T.activity
ifc.json_metadata.timesMid(end-20:end)
ifc
ifc.img = [ifc.img(1:470), T.activity'];
ifc.json_metadata.timesMid = [ifc.json_metadata.timesMid(1:470), T.t_elapsed'];
ifc.json_metadata.timesMid = [ifc.json_metadata.timesMid(1:470), T.t_elapsed];
size(ifc.json_metadata.timesMid)
ifc.json_metadata.timesMid = [ifc.json_metadata.timesMid(1:470); T.t_elapsed];
ifc.json_metadata.times = [ifc.json_metadata.timesMid(1:470); (T.t_elapsed - 1)];
ifc.json_metadata.taus(end-10:end)
ifc.json_metadata.taus = ifc.json_metadata.taus(1:474);
plot(ifc)
ic = mlfourd.ImagingContext2(ifc); plot(ic)
ifc
ifc.save
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108237/ses-20221031113804/pet')
ic = mlfourd.ImagingContext2("sub-108237_ses-20221031113804_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic)
+ifc = ic.imagingFormat;
ifc = ic.imagingFormat;
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,10,31))
t_drawn = rm.countsFdg.TIMEDRAWN_Hh_mm_ss;
activity = rm.countsFdg.DECAYCorrSpecificActivity_KBq_mL;
t_elapsed = seconds(t_drawn - datetime(2022,10,31,11,38,4, TimeZone="local")) + 17
T = table(t_drawn, t_elapsed, activity)
t_elapsed = seconds(t_drawn - datetime(2022,10,31,12,8,4, TimeZone="local")) + 17
T = table(t_drawn, t_elapsed, activity)
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,10,31))
open mlpet.CCIRRadMeasurements
rm.tracerAdmin
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,10,31))
rm.tracerAdmin
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,10,31))
rm.countsFdgVenous
t_drawn2 = rm.countsFdgVenous.TIMEDRAWN_Hh_mm_ss;
activity2 = rm.countsFdgVenous.DECAYCorrSpecificActivity_KBq_mL;
t_elapsed2 = seconds(t_drawn2 - datetime(2022,10,31,12,8,4, TimeZone="local")) + 17
open mlpet.CCIRRadMeasurements
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,10,31))
T
U = this.readtable(fqfn, 'Radiation Counts Log - Runs-2', 1, 1)
U = this.readtable(fqfn, 'Radiation Counts Log - Runs-2', 1, 1, 'exceldatenum')
U = this.readtable(fqfn, 'Radiation Counts Log - Runs-2', 1, 1, 'text')
U = this.readtable(fqfn, 'Radiation Counts Log - Runs-2', 1, 1, 'datetime')
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,10,31))
t_elapsed2 = seconds(t_drawn2 - datetime(2022,10,31,12,8,4, TimeZone="local")) + 17
t_drawn2 = rm.countsFdgVenous.TIMEDRAWN_Hh_mm_ss;
-t_elapsed2 = seconds(t_drawn2 - datetime(2022,10,31,12,8,4, TimeZone="local")) + 17
t_elapsed2 = seconds(t_drawn2 - datetime(2022,10,31,12,8,4, TimeZone="local")) + 17
t_elapsed2 = seconds(t_drawn2 - datetime(2022,10,31,12,8,4)) + 17
T = table(t_drawn, t_elapsed, activity)
T2 = table(t_drawn2, t_elapsed2, activity2)
U = [T; T2]
T2.Properties.VariableNames = T.Properties.VariableNames;
U = [T; T2]
T2.t_drawn.TimeZone = "local"
U = [T; T2]
U([2,4,6,8,10,12,14,end],:) = [];
U = sortrows(U, "t_elapsed")
ifc
ifc.json_metadata.timesMid(1:10)
ifc.json_metadata.timesMid(end)
ifc.img(332)
ifc
ifc.img = [ifc.img, asrow(U.activity(2:end))];
ifc.json_metadata.timesMid = [ifc.json_metadata.timesMid; U.t_elapsed(2:end)];
ifc.json_metadata.times(1:10)
ifc.json_metadata.times = [ifc.json_metadata.times; U.t_elapsed(2:end) - 1];
ifc
ifc.json_metadata
ifc.json_metadata.taus = [ifc.json_metadata.taus; ones(7,1)]
ifc.json_metadata
ifc.json_metadata.taus(end-10:end)
ifc.json_metadata.times(end-10:end)
ifc.json_metadata.timesMid(end-10:end)
ifc
ifc.save
ic = mlfourd.ImagingContext2(ifc); plot(ic)
figure; plot(ifc.img(end-20:end))
%-- 2/6/24, 15:56 --%
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108237/ses-20221031113804/pet')
ic = mlfourd.ImagingContext2("sub-108237_ses-20221031113804_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
ifc = ic.imagingFormat;
plot(ic)
ifc.json_metadata.timesMid(end-15:end)
ifc.img(406:end) = ifc.img(406:end) .* 2^(ifc.json_metadata.timesMid(406:end)/6586.272)
ifc.img(406:end) = ifc.img(406:end) .* 2^(ifc.json_metadata.timesMid(406:end)'/6586.272)
ifc.img(406:end) = ifc.img(406:end) .* 2^(asrow(ifc.json_metadata.timesMid(406:end))/6586.272)
ifc.img(406:end) = ifc.img(406:end) .* 2.^(asrow(ifc.json_metadata.timesMid(406:end))/6586.272)
ic = mlfourd.ImagingContext2(ifc); plot(ic)
ifc.img(406:end) = 1e3*ifc.img(406:end);
ic = mlfourd.ImagingContext2(ifc); plot(ic)
ic = mlfourd.ImagingContext2("sub-108237_ses-20221031113804_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
ifc = ic.imagingFormat;
plot(ic)
figure; plot(ifc.img(400:end))
ifc.json_metadata.timesMid(end-20:end)
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,10,31))
t_drawn2 = rm.countsFdgVenous.TIMEDRAWN_Hh_mm_ss;
activity2 = rm.countsFdgVenous.DECAYCorrSpecificActivity_KBq_mL;
t_elapsed2 = seconds(t_drawn2 - datetime(2022,10,31,12,8,4)) + 17
T2 = table(t_drawn2, t_elapsed2, activity2)
t_drawn = rm.countsFdg.TIMEDRAWN_Hh_mm_ss;
activity = rm.countsFdg.DECAYCorrSpecificActivity_KBq_mL;
t_elapsed = seconds(t_drawn - datetime(2022,10,31,12,8,4)) + 17
t_elapsed = seconds(t_drawn - datetime(2022,10,31,12,8,4), TimeZone="local") + 17
t_elapsed = seconds(t_drawn - datetime(2022,10,31,12,8,4, TimeZone="local")) + 17
T = table(t_drawn, t_elapsed, activity)
T2
T2.Properties.VariableNames = T.Properties.VariableNames;
U = sortrows(U, "t_elapsed")
U = [T; T2]
T.t_drawn.TimeZone = "local";
U = [T; T2]
T2.t_drawn.TimeZone = "local";
U = [T; T2]
U = sortrows(U, "t_elapsed")
U([2,4,6,8,10,12,14,end],:) = [];
U
U.activity = 1e3*U.activity .* 2.^(U.t_elapsed/6586.272)
ifc.json_metadata.timesMid(end-9:end)
ifc.img(end-9:end)
ifc.img(end-30:end)
ifc
ifc.json_metadata
ifc.img(407:end) = asrow(U.activity(2:end));
ic = mlfourd.ImagingContext2(ifc); plot(ic)
ifc.img(401:end)
ifc.img(396:end)
ifc.img(391:end)
ifc.img(393:406) = nan;
select = ~isnan(ifc.img);
ifc.img = ifc.img(select);
ifc.json_metadata.timesMid = ifc.json_metadata.timesMid(select);
ifc.json_metadata.times = ifc.json_metadata.times(select);
ifc.json_metadata.taus = ifc.json_metadata.taus(select);
ic = mlfourd.ImagingContext2(ifc); plot(ic)
ifc.save
close all
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108293/ses-20210421155709/pet')
ic = mlfourd.ImagingContext2("sub-108293_ses-20210421155709_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic)
ifc = ic.imagingFormat;
ifc.img(end-10:end)
ifc.json_metadata.timesMid(end-10:end)
ifc.img(end-3:end) = ifc.img(end-3:end) .* 2.^(asrow(ifc.json_metadata.timesMid(end-3:end)) / 6586.272)
ic = mlfourd.ImagingContext2(ifc); plot(ic)
ifc.save
close all
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2023,2,27))
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2022,12,7))
open mlpet.CCIRRadMeasurements
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2023,2,27))
rm.countsFdgVenous
rm = mlpet.CCIRRadMeasurements.createFromDate(datetime(2023,2,27))
cd('/Users/jjlee/Documents/CCIRRadMeasurements')
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 20230227.xlsx")
this.tracerAdmin
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 20230227.xlsx")
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 2022dec7.xlsx")
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 20230227.xlsx")
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 2022dec7.xlsx")
%-- 2/6/24, 19:33 --%
cd('/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn')
nB = mglob("NumBases*")
for idx = 1:20
ensuredir(fullfile(nB(idx), "components")); end
comp = mglob("NumBases*/components")
comp = natsort(comp)
cd(comp(1))
nmfr = mladni.NMFRegression(); nmfr.build_basis_argmax()
ic = ans;
for idx = 2:20
pwd0 = pushd(comp(idx)); nmfr = mladni.NMFRegression(); nmfr.build_basis_argmax(); popd(pwd0); end
mlsystem.Newcl.makecl("mladni", "NMFHierarchies", inheritance="handle")
ic
cd('/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn/NumBases2/OPNMF/niiImg')
ic = mlfourd.ImagingContext2("Basis_argmax_brain_mask.nii")
ic = ic.binarized
dipsum(ic)
91*109*91
num2str([1;2;3])
num2string([1;2;3])
double2str([1;2;3])
mat2str([1;2;3])
string(num2str([1;2;3]))
s = string(num2str([1;2;3]));
s(s == "1") = "found"
num2str(magic(3))
sortrows(magic(4), 4)
magic(4)
eps('single')
ic1 = mlfourd.ImagingContext2("Basis_1.nii");
ic2 = mlfourd.ImagingContext2("Basis_2.nii");
dipsum(ic1)
dipsum(ic2)
open mlfourd.MaskingTool
ic1.maskedMaths(ic1.numgt(1e-6), @mean)
ifc1 = ic1.imagingFormat;
ifc2 = ic2.imagingFormat;
dipmean(ifc1.img(ifc1.img > 1e-6))
dipmean(ifc2.img(ifc2.img > 1e-6))
nmfh = mladni.NMFHierarchies
nmfh.build_argmax_maps
nmfh = mladni.NMFHierarchies
nmfh.build_argmax_maps
nmfh = mladni.NMFHierarchies
nmfh.build_argmax_maps
nmfh = mladni.NMFHierarchies
nmfh.build_argmax_maps
nmfh = mladni.NMFHierarchies
nmfh.build_argmax_maps
nmfh = mladni.NMFHierarchies
nmfh.build_argmax_maps
T
load('/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn/NumBases24/components/NMFCovariates_table_covariates_1stscan_longitudinal.mat')
t
nmfr = mladni.NMFRadar
nmfr.table_patt_weighted_fdg
nmfr.table_patt_weighted_fdg.indices_bases
nmfr.table_patt_weighted_fdg.indices_bases'
nmfh = mladni.NMFHierarchies
nmfh.build_argmax_maps
nmfh = mladni.NMFHierarchies
nmfh.build_table_for_ggalluvial
nmfh = mladni.NMFHierarchies
nmfh.build_table_for_ggalluvial
"Span-"+asrow(string(num2str(this.selected_spans)))
asrow(string(num2str(this.selected_spans)))
string(num2str(this.selected_spans))
num2str(this.selected_spans)
string(this.selected_spans)
"Span-"+string(this.selected_spans)
nmfh = mladni.NMFHierarchies
nmfh.build_table_for_ggalluvial
"Span-"+asrow(string(this.selected_spans))
T
this.selected_spans
nmfh = mladni.NMFHierarchies
nmfh.build_table_for_ggalluvial
A(1:6,:)
imshow(A)
close
max(A)
dipmax(A)
dipmin(A)
[min(A(:,1) max(A(:,1))]
[min(A(:,1)) max(A(:,1))]
[min(A(:,2)) max(A(:,2))]
[min(A(:,end)) max(A(:,end))]
unique(tag)
tmp = string(A(:, end));
tmp(1:10)
this.mu_suvr_24
[max(suvr) min(suvr)\
[max(suvr) min(suvr)]
histogram(suvr)
close
nmfh = mladni.NMFHierarchies
nmfh.build_table_for_ggalluvial
T
unique(T.Tag)
nmfh = mladni.NMFHierarchies
nmfh.build_table_for_ggalluvial
open mlkinetic.QuadraticRaichle1983Model
open mlkinetics.QuadraticRaichle1983Model
open mlkinetics.QuadraticMintun1984Model
open mloxygen.QuadraticRaichle1983Model
open mloxygen.Raichle1983Model
cd('/Users/jjlee/Documents/CCIRRadMeasurements')
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 2022dec7.xlsx")
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108250/ses-20221207104909/pet')
ic = mlfourd.ImagingContext2("sub-108250_ses-20221207104909_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic)
ifc = ic.imagingFormat;
ifc.json_metadata.timesMid(end-15:end)
3110/3600
50/60
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 2022dec7.xlsx")
cd('/Users/jjlee/Documents/CCIRRadMeasurements')
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 2022dec7.xlsx")
ifc
2^(50/110)
5*2^(50/110)
2^(20/110)
6*2^(50/110)
5.7*2^(20/110)
60*3110/3600
close
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108306/ses-20230227115809/pet')
ic = mlfourd.ImagingContext2("sub-108306_ses-20230227115809_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic)
cd('/Users/jjlee/Documents/CCIRRadMeasurements')
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 20220227.xlsx")
rm = mlpet.CCIRRadMeasurements.createFromFilename("CCIRRadMeasurements 20230227.xlsx")
6*2^(25/110)
3600*3/6
ifc = ic.imagingFormat;
ifc.json_metadata.timesMid(end-15:end)
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108306/ses-20230227115809/pet')
t_drawn2 = rm.countsFdgVenous.TIMEDRAWN_Hh_mm_ss;
activity2 = rm.countsFdgVenous.DECAYCorrSpecificActivity_KBq_mL;
t_elapsed2 = seconds(t_drawn2 - datetime(2023,2,27,11,38,53)) + 32
T2 = table(t_drawn2, t_elapsed2, activity2)
t_drawn2 = rm.countsFdgVenous.TIMEDRAWN_Hh_mm_ss;
activity2 = rm.countsFdgVenous.DECAYCorrSpecificActivity_KBq_mL;
t_elapsed2 = seconds(t_drawn2 - datetime(2023,2,27,11,38,53))
T2 = table(t_drawn2, t_elapsed2, activity2)
T2.activity2 = 1e3*T2.activity2 .* 2.^(T2.t_elapsed2/6586.272)
ifc
figure; plot(ifc.img)
T2.activity2'
img = ifc.img; img = [img, T2.activity2'];
img(1695:end)
size(img)
img(1595:end)
j = ifc.json_metadata
j.taus(end)
j.times(end-10:end)
j.timesMid(end-10:end)
j.timesMid = [j.timesMid; T2.t_elapsed2];
j.times = [j.times; T2.t_elapsed2 - 1];
j
j.taus = [j.taus; ones(6,1)];
j
j.taus(end-10:end)
ifc.json_metadata = j;
ifc
ifc.img = img
ifc.save
ic = mlfourd.ImagingContext2(ifc); plot(ic)
close all
ic = mlfourd.ImagingContext2("sub-108306_ses-20230227115809_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.nii.gz")
plot(ic)
t2 = [1,2,3,4,5]; t1 = [1,2,3];
[a,b] = t2 > t1(end)
[a,b] = max(t2 > t1(end))
ic1 = mlfourd.ImagingContext2("sub-108306_ses-20230227115809_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.nii.gz")
ic2 = mlfourd.ImagingContext2("sub-108306_ses-20230227115809_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic1)
plot(ic2)
ic = mlvg.Lee2024.append_activity_densities(ic1, ic2)
ic1
ic
close all
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108250/ses-20221207104909/pet')
ic1 = mlfourd.ImagingContext2("sub-108250_ses-20221207104909_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.nii.gz")
ic2 = mlfourd.ImagingContext2("sub-108250_ses-20221207104909_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic1); plot(ic2);
ic = mlvg.Lee2024.append_activity_densities(ic1, ic2)
close all
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108237/ses-20221031113804/pet')
ic1 = mlfourd.ImagingContext2("sub-108237_ses-20221031113804_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.nii.gz")
ic2 = mlfourd.ImagingContext2("sub-108237_ses-20221031113804_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic1); plot(ic2);
ic = mlvg.Lee2024.append_activity_densities(ic1, ic2)
close all
cd('/Volumes/PrecunealSSD/Singularity/CCIR_01211/derivatives/sub-108293/ses-20210421155709/pet')
ic1 = mlfourd.ImagingContext2("sub-108293_ses-20210421155709_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc_dynesty-RadialArtery-ideal.nii.gz")
ic2 = mlfourd.ImagingContext2("sub-108293_ses-20210421155709_trc-fdg_proc-TwiliteKit-do-make-input-func-nomodel_inputfunc.nii.gz")
plot(ic1); plot(ic2);
ic = mlvg.Lee2024.append_activity_densities(ic1, ic2)
close all