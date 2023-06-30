module load dSQ
#!/bin/bash
numSeeds=1
calA0=(1.0)
phi_arr=(0.74 0.63 0.53)
att_arr=(0 0.001 0.01 0.1)
v0_arr=(0.005)
kecm_arr=(0.006 0.06 0.6)
koff_arr=(1.0 10.0 100.0)
rm joblist_psm_att_v0_koff.txt
for phi in ${phi_arr[@]}; do
  for att in ${att_arr[@]}; do 
    for v0 in ${v0_arr[@]}; do
      for koff in ${koff_arr[@]}; do
          for kecm in ${kecm_arr[@]}; do
              echo bash bash/cells/submit_psm.sh 40 20 $calA0 $phi $att 0.0 $v0 1.0 $kecm $koff 1000 pi_ohern,day 0-12:00:00 $numSeeds 1 >> joblist_psm_att_v0_koff.txt
          done
      done
    done
  done
done

dsq --job-file joblist_psm_att_v0_koff.txt   --mem-per-cpu 4g -t 1:00:00 --mail-type NONE --submit --partition scavenge --suppress-stats-file  -o /dev/null

bash bash/cells/submit_psm.sh 40 20 1.05 0.9 0 0.0 0.05 50.0 1.0 10.0 1000 pi_ohern,day 0-4:00:00 1 1

rsync -rav --inplace --progress at965@transfer-grace.hpc.yale.edu:/gpfs/gibbs/pi/ohern/at965/dpm/psm /mnt/c/Users/atata/projects/dpm/pipeline/cells/. 

close all; clear;
calA0_arr = ["1.0"];
att_arr = ["0" "0.001" "0.01" "0.1"];
phi_arr = ["0.74" "0.63" "0.53"]
v0_arr = ["0.005"];
k_ecm_arr = ["0.006" "0.06" "0.6"];
k_off_arr = ["1.0" "10.0" "100.0"];
for ii=1:length(calA0_arr)
    for jj=1:length(phi_arr)
        for kk=1:length(att_arr)
            for ll=1:length(v0_arr)
                for mm=1:length(k_ecm_arr)
                  for nn=1:length(k_off_arr)
                    drawCellSim("40", calA0_arr(ii), phi_arr(jj), att_arr(kk), v0_arr(ll), k_ecm_arr(mm), k_off_arr(nn))
                  end
                end
            end
        end
    end
end

MATLAB Log File: C:\Users\atata\AppData\Local\Temp\matlab_crash_dump.17028-1

------------------------------------------------
MATLAB Log File
------------------------------------------------ 


--------------------------------------------------------------------------------
             Access violation detected at 2023-06-30 09:40:22 -0400
--------------------------------------------------------------------------------

Configuration:
  Crash Decoding           : Disabled - No sandbox or build area path
  Crash Mode               : continue (default)
  Default Encoding         : UTF-8
  Deployed                 : false
  Graphics Driver          : NVIDIA Corporation NVIDIA GeForce RTX 3060 Ti/PCIe/SSE2 Version 4.6.0 NVIDIA 516.40
  Graphics card 1          : NVIDIA ( 0x10de ) NVIDIA GeForce RTX 3060 Ti Version 31.0.15.1640 (2022-6-6)
  Java Version             : Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
  MATLAB Architecture      : win64
  MATLAB Entitlement ID    : 873034
  MATLAB Root              : C:\Program Files\MATLAB\R2022a
  MATLAB Version           : 9.12.0.1927505 (R2022a) Update 1
  OpenGL                   : hardware
  Operating System         : Microsoft Windows 11 Home
  Process ID               : 17028
  Processor ID             : x86 Family 23 Model 113 Stepping 0, AuthenticAMD
  Session Key              : 4c97b570-782b-49f1-ac37-8ec8e5126b59
  Window System            : Version 10.0 (Build 22621)

Fault Count: 1


Abnormal termination:
Access violation

Current Thread: 'MCR 0 interpreter thread' id 9480

Register State (from fault):
  RAX = 000000eb2cdfa9e8  RBX = 000000000015f8c4
  RCX = 0000002000000000  RDX = 000002c5938ab060
  RSP = 000000eb2cdfa970  RBP = 000002c5938ab060
  RSI = 0000000002719c40  RDI = 000002c5938ab060
 
   R8 = 000002c52f711720   R9 = 0000000000000008
  R10 = 0000000000000001  R11 = 000000eb2cdfa8d8
  R12 = 000000eb2cdfafc0  R13 = 0000000000000001
  R14 = 00007ffceda66230  R15 = 0000000000000001
 
  RIP = 00007ffceda6623a  EFL = 00010202
 
   CS = 0033   FS = 0053   GS = 002b

Stack Trace (from fault):
[  0] 0x00007ffceda6623a                            bin\win64\pgo\libmx.dll+00025146 matrix::detail::noninlined::mx_array_api::mxDestroyArray+00000010
[  1] 0x00007ffceda6445f                            bin\win64\pgo\libmx.dll+00017503 matrix::detail::noninlined::mx_array_api::mxGetFieldNumber+00000767
[  2] 0x00007ffceda74c57                            bin\win64\pgo\libmx.dll+00085079 matrix::detail::noninlined::mx_array_api::mxMallocExCheck+00002055
[  3] 0x00007ffceda643d7                            bin\win64\pgo\libmx.dll+00017367 matrix::detail::noninlined::mx_array_api::mxGetFieldNumber+00000631
[  4] 0x00007ffceda642ee                            bin\win64\pgo\libmx.dll+00017134 matrix::detail::noninlined::mx_array_api::mxGetFieldNumber+00000398
[  5] 0x00007ffceda74fd4                            bin\win64\pgo\libmx.dll+00085972 matrix::detail::noninlined::mx_array_api::mxSubscriptedDelete+00000404
[  6] 0x00007ffce666c0e0                 bin\win64\pgo\libmwlxeindexing.dll+00049376 MathWorks::lxe::paren_delete_pointer_cell+00000320
[  7] 0x00007ffce75d79b8 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+04815288 MathWorks::lxe::printLxeProfStatsForFeature+00279032
[  8] 0x00007ffce75f7449 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+04944969 MathWorks::lxe::printLxeProfStatsForFeature+00408713
[  9] 0x00007ffce73ebbab C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+02800555 MathWorks::lxe::LXEConstants::IsIfElse+00522971
[ 10] 0x00007ffce73f0074 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+02818164 MathWorks::lxe::LXEConstants::IsY+00012660
[ 11] 0x00007ffce73ed601 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+02807297 MathWorks::lxe::LXEConstants::IsY+00001793
[ 12] 0x00007ffce73f11b5 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+02822581 MathWorks::lxe::LXEConstants::IsY+00017077
[ 13] 0x00007ffce73f1463 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+02823267 MathWorks::lxe::LXEConstants::IsY+00017763
[ 14] 0x00007ffce73f0db8 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+02821560 MathWorks::lxe::LXEConstants::IsY+00016056
[ 15] 0x00007ffce73ecda6 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+02805158 MathWorks::lxe::LxeTypes::GetTypeXvalueOf+00000598
[ 16] 0x00007ffce74bb3e5 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+03650533 MathWorks::lxe::SetEngineImplUnlinkFlagForAllCallsOnStack+00043317
[ 17] 0x00007ffce74bef60 C:\Program Files\MATLAB\R2022a\bin\win64\m_lxe.dll+03665760 MathWorks::lxe::SetEngineImplUnlinkFlagForAllCallsOnStack+00058544
[ 18] 0x00007ffcea4df830 C:\Program Files\MATLAB\R2022a\bin\win64\libmwlxemainservices.dll+01177648 FeatureTestObservableWorkspace+00186352
[ 19] 0x00007ffcea44bc7e C:\Program Files\MATLAB\R2022a\bin\win64\libmwlxemainservices.dll+00572542 mwboost::archive::detail::pointer_oserializer<mwboost::archive::xml_oarchive,MathWorks::lxe::PreLineExecutionEvent>::save_object_ptr+00069966
[ 20] 0x00007ffcea44d5bf C:\Program Files\MATLAB\R2022a\bin\win64\libmwlxemainservices.dll+00579007 mwboost::archive::detail::pointer_oserializer<mwboost::archive::xml_oarchive,MathWorks::lxe::PreLineExecutionEvent>::save_object_ptr+00076431
[ 21] 0x00007ffcea4a8cf6 C:\Program Files\MATLAB\R2022a\bin\win64\libmwlxemainservices.dll+00953590 mwboost::archive::detail::pointer_oserializer<mwboost::archive::xml_oarchive,MathWorks::lxe::PreLineExecutionEvent>::save_object_ptr+00451014
[ 22] 0x00007ffcea4aaca7 C:\Program Files\MATLAB\R2022a\bin\win64\libmwlxemainservices.dll+00961703 mwboost::archive::detail::pointer_oserializer<mwboost::archive::xml_oarchive,MathWorks::lxe::PreLineExecutionEvent>::save_object_ptr+00459127
[ 23] 0x00007ffcea4aa050 C:\Program Files\MATLAB\R2022a\bin\win64\libmwlxemainservices.dll+00958544 mwboost::archive::detail::pointer_oserializer<mwboost::archive::xml_oarchive,MathWorks::lxe::PreLineExecutionEvent>::save_object_ptr+00455968
[ 24] 0x00007ffceaa725b9 C:\Program Files\MATLAB\R2022a\bin\win64\libmwbridge.dll+00140729 mnGetPrompt+00025257
[ 25] 0x00007ffceab275c6   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00488902 iqm::Iqm::instance+00001190
[ 26] 0x00007ffceab76a60   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00813664 iqm::UserEvalPlugin::execute+00001360
[ 27] 0x00007ffc8c6f2b75 C:\Program Files\MATLAB\R2022a\bin\win64\nativejmi.dll+00338805 nativejmi::JmiIIP::serializeExplicit+00008005
[ 28] 0x00007ffceab49882   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00628866 iqm::Iqm::setupIqmFcnPtrs+00094722
[ 29] 0x00007ffceab533ca   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00668618 iqm::Iqm::setupIqmFcnPtrs+00134474
[ 30] 0x00007ffceab1af4b   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00438091 iqm::Iqm::create+00009067
[ 31] 0x00007ffceaa62de5 C:\Program Files\MATLAB\R2022a\bin\win64\libmwbridge.dll+00077285 ioReadLine+00000501
[ 32] 0x00007ffceaa62bb5 C:\Program Files\MATLAB\R2022a\bin\win64\libmwbridge.dll+00076725 ioReadLine+00000165
[ 33] 0x00007ffceaa72970 C:\Program Files\MATLAB\R2022a\bin\win64\libmwbridge.dll+00141680 mnGetCommandLineBuffer+00000288
[ 34] 0x00007ffceaa72e02 C:\Program Files\MATLAB\R2022a\bin\win64\libmwbridge.dll+00142850 mnParser+00000466
[ 35] 0x00007ffceac87226   C:\Program Files\MATLAB\R2022a\bin\win64\mcr.dll+00422438 mcr_set_enableReadingFromStdin+00013622
[ 36] 0x00007ffceac389a3   C:\Program Files\MATLAB\R2022a\bin\win64\mcr.dll+00100771 mcrFunctionSignature::set_signature+00079731
[ 37] 0x00007ffceac56230   C:\Program Files\MATLAB\R2022a\bin\win64\mcr.dll+00221744 mwboost::archive::codecvt_null<wchar_t>::`default constructor closure'+00017584
[ 38] 0x00007ffceab71e0a   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00794122 iqm::PackagedTaskPlugin::execute+00000074
[ 39] 0x00007ffceac819ad   C:\Program Files\MATLAB\R2022a\bin\win64\mcr.dll+00399789 services::lmgr::exception::LicensingStartupException::~LicensingStartupException+00006861
[ 40] 0x00007ffceab49882   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00628866 iqm::Iqm::setupIqmFcnPtrs+00094722
[ 41] 0x00007ffceab1bd86   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00441734 iqm::Iqm::create+00012710
[ 42] 0x00007ffceab1b472   C:\Program Files\MATLAB\R2022a\bin\win64\iqm.dll+00439410 iqm::Iqm::create+00010386
[ 43] 0x00007ffceac7072c   C:\Program Files\MATLAB\R2022a\bin\win64\mcr.dll+00329516 mcrInstantiationError::operator=+00010380
[ 44] 0x00007ffceac71165   C:\Program Files\MATLAB\R2022a\bin\win64\mcr.dll+00332133 mcrInstantiationError::operator=+00012997
[ 45] 0x00007ffceac6e930   C:\Program Files\MATLAB\R2022a\bin\win64\mcr.dll+00321840 mcrInstantiationError::operator=+00002704
[ 46] 0x00007ffd06d285da C:\Program Files\MATLAB\R2022a\bin\win64\mwboost_thread-vc142-mt-x64-1_75.dll+00034266 mwboost::thread::swap+00000074
[ 47] 0x00007ffd2a0c9363                   C:\WINDOWS\System32\ucrtbase.dll+00168803 recalloc+00000163
[ 48] 0x00007ffd2b9126ad                   C:\WINDOWS\System32\KERNEL32.DLL+00075437 BaseThreadInitThunk+00000029
[ 49] 0x00007ffd2ca6a9f8                      C:\WINDOWS\SYSTEM32\ntdll.dll+00371192 RtlUserThreadStart+00000040

