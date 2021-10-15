//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 13 21:27:30 2021 by ROOT version 6.14/04
// from TTree T/Hall A Analyzer Output DST
// found on file: replayed_simdigtest_2_20211004.root
//////////////////////////////////////////////////////////

#ifndef gmn_rec_tree_h
#define gmn_rec_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class gmn_rec_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Ndata_MC_bbgemhit_edep;
   Double_t        MC_bbgemhit_edep[68];   //[Ndata.MC.bbgemhit_edep]
   Int_t           Ndata_MC_bbgemhit_mid;
   Double_t        MC_bbgemhit_mid[68];   //[Ndata.MC.bbgemhit_mid]
   Int_t           Ndata_MC_bbgemhit_pid;
   Double_t        MC_bbgemhit_pid[68];   //[Ndata.MC.bbgemhit_pid]
   Int_t           Ndata_MC_bbgemhit_plane;
   Double_t        MC_bbgemhit_plane[68];   //[Ndata.MC.bbgemhit_plane]
   Int_t           Ndata_MC_bbgemhit_tid;
   Double_t        MC_bbgemhit_tid[68];   //[Ndata.MC.bbgemhit_tid]
   Int_t           Ndata_MC_bbgemhit_x;
   Double_t        MC_bbgemhit_x[68];   //[Ndata.MC.bbgemhit_x]
   Int_t           Ndata_MC_bbgemhit_y;
   Double_t        MC_bbgemhit_y[68];   //[Ndata.MC.bbgemhit_y]
   Int_t           Ndata_MC_bbtrack_dx;
   Double_t        MC_bbtrack_dx[11];   //[Ndata.MC.bbtrack_dx]
   Int_t           Ndata_MC_bbtrack_dy;
   Double_t        MC_bbtrack_dy[11];   //[Ndata.MC.bbtrack_dy]
   Int_t           Ndata_MC_bbtrack_mid;
   Double_t        MC_bbtrack_mid[11];   //[Ndata.MC.bbtrack_mid]
   Int_t           Ndata_MC_bbtrack_nhits;
   Double_t        MC_bbtrack_nhits[11];   //[Ndata.MC.bbtrack_nhits]
   Int_t           Ndata_MC_bbtrack_p;
   Double_t        MC_bbtrack_p[11];   //[Ndata.MC.bbtrack_p]
   Int_t           Ndata_MC_bbtrack_pid;
   Double_t        MC_bbtrack_pid[11];   //[Ndata.MC.bbtrack_pid]
   Int_t           Ndata_MC_bbtrack_tid;
   Double_t        MC_bbtrack_tid[11];   //[Ndata.MC.bbtrack_tid]
   Int_t           Ndata_MC_bbtrack_x;
   Double_t        MC_bbtrack_x[11];   //[Ndata.MC.bbtrack_x]
   Int_t           Ndata_MC_bbtrack_y;
   Double_t        MC_bbtrack_y[11];   //[Ndata.MC.bbtrack_y]
   Int_t           Ndata_MC_pt_clustsz;
   Double_t        MC_pt_clustsz[1];   //[Ndata.MC.pt.clustsz]
   Int_t           Ndata_MC_pt_deflect;
   Double_t        MC_pt_deflect[1];   //[Ndata.MC.pt.deflect]
   Int_t           Ndata_MC_pt_deltaE;
   Double_t        MC_pt_deltaE[1];   //[Ndata.MC.pt.deltaE]
   Int_t           Ndata_MC_pt_hitres;
   Double_t        MC_pt_hitres[1];   //[Ndata.MC.pt.hitres]
   Int_t           Ndata_MC_pt_nfound;
   Double_t        MC_pt_nfound[1];   //[Ndata.MC.pt.nfound]
   Int_t           Ndata_MC_pt_p;
   Double_t        MC_pt_p[1];   //[Ndata.MC.pt.p]
   Int_t           Ndata_MC_pt_ph;
   Double_t        MC_pt_ph[1];   //[Ndata.MC.pt.ph]
   Int_t           Ndata_MC_pt_phdir;
   Double_t        MC_pt_phdir[1];   //[Ndata.MC.pt.phdir]
   Int_t           Ndata_MC_pt_phi;
   Double_t        MC_pt_phi[1];   //[Ndata.MC.pt.phi]
   Int_t           Ndata_MC_pt_plane;
   Double_t        MC_pt_plane[1];   //[Ndata.MC.pt.plane]
   Int_t           Ndata_MC_pt_r;
   Double_t        MC_pt_r[1];   //[Ndata.MC.pt.r]
   Int_t           Ndata_MC_pt_status;
   Double_t        MC_pt_status[1];   //[Ndata.MC.pt.status]
   Int_t           Ndata_MC_pt_th;
   Double_t        MC_pt_th[1];   //[Ndata.MC.pt.th]
   Int_t           Ndata_MC_pt_thdir;
   Double_t        MC_pt_thdir[1];   //[Ndata.MC.pt.thdir]
   Int_t           Ndata_MC_pt_theta;
   Double_t        MC_pt_theta[1];   //[Ndata.MC.pt.theta]
   Int_t           Ndata_MC_pt_time;
   Double_t        MC_pt_time[1];   //[Ndata.MC.pt.time]
   Int_t           Ndata_MC_pt_tof;
   Double_t        MC_pt_tof[1];   //[Ndata.MC.pt.tof]
   Int_t           Ndata_MC_pt_trkres;
   Double_t        MC_pt_trkres[1];   //[Ndata.MC.pt.trkres]
   Int_t           Ndata_MC_pt_type;
   Double_t        MC_pt_type[1];   //[Ndata.MC.pt.type]
   Int_t           Ndata_MC_pt_x;
   Double_t        MC_pt_x[1];   //[Ndata.MC.pt.x]
   Int_t           Ndata_MC_pt_y;
   Double_t        MC_pt_y[1];   //[Ndata.MC.pt.y]
   Int_t           Ndata_bb_gem_hit_ADCU;
   Double_t        bb_gem_hit_ADCU[12];   //[Ndata.bb.gem.hit.ADCU]
   Int_t           Ndata_bb_gem_hit_ADCV;
   Double_t        bb_gem_hit_ADCV[12];   //[Ndata.bb.gem.hit.ADCV]
   Int_t           Ndata_bb_gem_hit_ADCasym;
   Double_t        bb_gem_hit_ADCasym[12];   //[Ndata.bb.gem.hit.ADCasym]
   Int_t           Ndata_bb_gem_hit_ADCavg;
   Double_t        bb_gem_hit_ADCavg[12];   //[Ndata.bb.gem.hit.ADCavg]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampU;
   Double_t        bb_gem_hit_ADCmaxsampU[12];   //[Ndata.bb.gem.hit.ADCmaxsampU]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampUclust;
   Double_t        bb_gem_hit_ADCmaxsampUclust[12];   //[Ndata.bb.gem.hit.ADCmaxsampUclust]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampV;
   Double_t        bb_gem_hit_ADCmaxsampV[12];   //[Ndata.bb.gem.hit.ADCmaxsampV]
   Int_t           Ndata_bb_gem_hit_ADCmaxsampVclust;
   Double_t        bb_gem_hit_ADCmaxsampVclust[12];   //[Ndata.bb.gem.hit.ADCmaxsampVclust]
   Int_t           Ndata_bb_gem_hit_ADCmaxstripU;
   Double_t        bb_gem_hit_ADCmaxstripU[12];   //[Ndata.bb.gem.hit.ADCmaxstripU]
   Int_t           Ndata_bb_gem_hit_ADCmaxstripV;
   Double_t        bb_gem_hit_ADCmaxstripV[12];   //[Ndata.bb.gem.hit.ADCmaxstripV]
   Int_t           Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U;
   Double_t        bb_gem_hit_BUILD_ALL_SAMPLES_U[12];   //[Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_U]
   Int_t           Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V;
   Double_t        bb_gem_hit_BUILD_ALL_SAMPLES_V[12];   //[Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_V]
   Int_t           Ndata_bb_gem_hit_CM_OR_U;
   Double_t        bb_gem_hit_CM_OR_U[12];   //[Ndata.bb.gem.hit.CM_OR_U]
   Int_t           Ndata_bb_gem_hit_CM_OR_V;
   Double_t        bb_gem_hit_CM_OR_V[12];   //[Ndata.bb.gem.hit.CM_OR_V]
   Int_t           Ndata_bb_gem_hit_ENABLE_CM_U;
   Double_t        bb_gem_hit_ENABLE_CM_U[12];   //[Ndata.bb.gem.hit.ENABLE_CM_U]
   Int_t           Ndata_bb_gem_hit_ENABLE_CM_V;
   Double_t        bb_gem_hit_ENABLE_CM_V[12];   //[Ndata.bb.gem.hit.ENABLE_CM_V]
   Int_t           Ndata_bb_gem_hit_Tavg;
   Double_t        bb_gem_hit_Tavg[12];   //[Ndata.bb.gem.hit.Tavg]
   Int_t           Ndata_bb_gem_hit_Utime;
   Double_t        bb_gem_hit_Utime[12];   //[Ndata.bb.gem.hit.Utime]
   Int_t           Ndata_bb_gem_hit_UtimeMaxStrip;
   Double_t        bb_gem_hit_UtimeMaxStrip[12];   //[Ndata.bb.gem.hit.UtimeMaxStrip]
   Int_t           Ndata_bb_gem_hit_Vtime;
   Double_t        bb_gem_hit_Vtime[12];   //[Ndata.bb.gem.hit.Vtime]
   Int_t           Ndata_bb_gem_hit_VtimeMaxStrip;
   Double_t        bb_gem_hit_VtimeMaxStrip[12];   //[Ndata.bb.gem.hit.VtimeMaxStrip]
   Int_t           Ndata_bb_gem_hit_ccor_clust;
   Double_t        bb_gem_hit_ccor_clust[12];   //[Ndata.bb.gem.hit.ccor_clust]
   Int_t           Ndata_bb_gem_hit_ccor_strip;
   Double_t        bb_gem_hit_ccor_strip[12];   //[Ndata.bb.gem.hit.ccor_strip]
   Int_t           Ndata_bb_gem_hit_deltat;
   Double_t        bb_gem_hit_deltat[12];   //[Ndata.bb.gem.hit.deltat]
   Int_t           Ndata_bb_gem_hit_eresidu;
   Double_t        bb_gem_hit_eresidu[12];   //[Ndata.bb.gem.hit.eresidu]
   Int_t           Ndata_bb_gem_hit_eresidv;
   Double_t        bb_gem_hit_eresidv[12];   //[Ndata.bb.gem.hit.eresidv]
   Int_t           Ndata_bb_gem_hit_isampmaxUclust;
   Double_t        bb_gem_hit_isampmaxUclust[12];   //[Ndata.bb.gem.hit.isampmaxUclust]
   Int_t           Ndata_bb_gem_hit_isampmaxUstrip;
   Double_t        bb_gem_hit_isampmaxUstrip[12];   //[Ndata.bb.gem.hit.isampmaxUstrip]
   Int_t           Ndata_bb_gem_hit_isampmaxVclust;
   Double_t        bb_gem_hit_isampmaxVclust[12];   //[Ndata.bb.gem.hit.isampmaxVclust]
   Int_t           Ndata_bb_gem_hit_isampmaxVstrip;
   Double_t        bb_gem_hit_isampmaxVstrip[12];   //[Ndata.bb.gem.hit.isampmaxVstrip]
   Int_t           Ndata_bb_gem_hit_layer;
   Double_t        bb_gem_hit_layer[12];   //[Ndata.bb.gem.hit.layer]
   Int_t           Ndata_bb_gem_hit_module;
   Double_t        bb_gem_hit_module[12];   //[Ndata.bb.gem.hit.module]
   Int_t           Ndata_bb_gem_hit_nstripu;
   Double_t        bb_gem_hit_nstripu[12];   //[Ndata.bb.gem.hit.nstripu]
   Int_t           Ndata_bb_gem_hit_nstripv;
   Double_t        bb_gem_hit_nstripv[12];   //[Ndata.bb.gem.hit.nstripv]
   Int_t           Ndata_bb_gem_hit_residu;
   Double_t        bb_gem_hit_residu[12];   //[Ndata.bb.gem.hit.residu]
   Int_t           Ndata_bb_gem_hit_residv;
   Double_t        bb_gem_hit_residv[12];   //[Ndata.bb.gem.hit.residv]
   Int_t           Ndata_bb_gem_hit_trackindex;
   Double_t        bb_gem_hit_trackindex[12];   //[Ndata.bb.gem.hit.trackindex]
   Int_t           Ndata_bb_gem_hit_u;
   Double_t        bb_gem_hit_u[12];   //[Ndata.bb.gem.hit.u]
   Int_t           Ndata_bb_gem_hit_umoment;
   Double_t        bb_gem_hit_umoment[12];   //[Ndata.bb.gem.hit.umoment]
   Int_t           Ndata_bb_gem_hit_usigma;
   Double_t        bb_gem_hit_usigma[12];   //[Ndata.bb.gem.hit.usigma]
   Int_t           Ndata_bb_gem_hit_ustriphi;
   Double_t        bb_gem_hit_ustriphi[12];   //[Ndata.bb.gem.hit.ustriphi]
   Int_t           Ndata_bb_gem_hit_ustriplo;
   Double_t        bb_gem_hit_ustriplo[12];   //[Ndata.bb.gem.hit.ustriplo]
   Int_t           Ndata_bb_gem_hit_ustripmax;
   Double_t        bb_gem_hit_ustripmax[12];   //[Ndata.bb.gem.hit.ustripmax]
   Int_t           Ndata_bb_gem_hit_v;
   Double_t        bb_gem_hit_v[12];   //[Ndata.bb.gem.hit.v]
   Int_t           Ndata_bb_gem_hit_vmoment;
   Double_t        bb_gem_hit_vmoment[12];   //[Ndata.bb.gem.hit.vmoment]
   Int_t           Ndata_bb_gem_hit_vsigma;
   Double_t        bb_gem_hit_vsigma[12];   //[Ndata.bb.gem.hit.vsigma]
   Int_t           Ndata_bb_gem_hit_vstriphi;
   Double_t        bb_gem_hit_vstriphi[12];   //[Ndata.bb.gem.hit.vstriphi]
   Int_t           Ndata_bb_gem_hit_vstriplo;
   Double_t        bb_gem_hit_vstriplo[12];   //[Ndata.bb.gem.hit.vstriplo]
   Int_t           Ndata_bb_gem_hit_vstripmax;
   Double_t        bb_gem_hit_vstripmax[12];   //[Ndata.bb.gem.hit.vstripmax]
   Int_t           Ndata_bb_gem_hit_xglobal;
   Double_t        bb_gem_hit_xglobal[12];   //[Ndata.bb.gem.hit.xglobal]
   Int_t           Ndata_bb_gem_hit_xlocal;
   Double_t        bb_gem_hit_xlocal[12];   //[Ndata.bb.gem.hit.xlocal]
   Int_t           Ndata_bb_gem_hit_yglobal;
   Double_t        bb_gem_hit_yglobal[12];   //[Ndata.bb.gem.hit.yglobal]
   Int_t           Ndata_bb_gem_hit_ylocal;
   Double_t        bb_gem_hit_ylocal[12];   //[Ndata.bb.gem.hit.ylocal]
   Int_t           Ndata_bb_gem_hit_zglobal;
   Double_t        bb_gem_hit_zglobal[12];   //[Ndata.bb.gem.hit.zglobal]
   Int_t           Ndata_bb_gem_n2Dhit_layer;
   Double_t        bb_gem_n2Dhit_layer[5];   //[Ndata.bb.gem.n2Dhit_layer]
   Int_t           Ndata_bb_gem_nclustu_layer;
   Double_t        bb_gem_nclustu_layer[5];   //[Ndata.bb.gem.nclustu_layer]
   Int_t           Ndata_bb_gem_nclustv_layer;
   Double_t        bb_gem_nclustv_layer[5];   //[Ndata.bb.gem.nclustv_layer]
   Int_t           Ndata_bb_gem_nstripsu_layer;
   Double_t        bb_gem_nstripsu_layer[5];   //[Ndata.bb.gem.nstripsu_layer]
   Int_t           Ndata_bb_gem_nstripsv_layer;
   Double_t        bb_gem_nstripsv_layer[5];   //[Ndata.bb.gem.nstripsv_layer]
   Int_t           Ndata_bb_gem_track_chi2ndf;
   Double_t        bb_gem_track_chi2ndf[3];   //[Ndata.bb.gem.track.chi2ndf]
   Int_t           Ndata_bb_gem_track_nhits;
   Double_t        bb_gem_track_nhits[3];   //[Ndata.bb.gem.track.nhits]
   Int_t           Ndata_bb_gem_track_x;
   Double_t        bb_gem_track_x[3];   //[Ndata.bb.gem.track.x]
   Int_t           Ndata_bb_gem_track_xp;
   Double_t        bb_gem_track_xp[3];   //[Ndata.bb.gem.track.xp]
   Int_t           Ndata_bb_gem_track_y;
   Double_t        bb_gem_track_y[3];   //[Ndata.bb.gem.track.y]
   Int_t           Ndata_bb_gem_track_yp;
   Double_t        bb_gem_track_yp[3];   //[Ndata.bb.gem.track.yp]
   Int_t           Ndata_bb_grinch_clus_adc;
   Double_t        bb_grinch_clus_adc[14];   //[Ndata.bb.grinch.clus.adc]
   Int_t           Ndata_bb_grinch_clus_size;
   Double_t        bb_grinch_clus_size[14];   //[Ndata.bb.grinch.clus.size]
   Int_t           Ndata_bb_grinch_clus_tf_mean;
   Double_t        bb_grinch_clus_tf_mean[14];   //[Ndata.bb.grinch.clus.tf_mean]
   Int_t           Ndata_bb_grinch_clus_tr_mean;
   Double_t        bb_grinch_clus_tr_mean[14];   //[Ndata.bb.grinch.clus.tr_mean]
   Int_t           Ndata_bb_grinch_clus_x_mean;
   Double_t        bb_grinch_clus_x_mean[14];   //[Ndata.bb.grinch.clus.x_mean]
   Int_t           Ndata_bb_grinch_clus_y_mean;
   Double_t        bb_grinch_clus_y_mean[14];   //[Ndata.bb.grinch.clus.y_mean]
   Int_t           Ndata_bb_grinch_hit_adc;
   Double_t        bb_grinch_hit_adc[170];   //[Ndata.bb.grinch.hit.adc]
   Int_t           Ndata_bb_grinch_hit_col;
   Double_t        bb_grinch_hit_col[170];   //[Ndata.bb.grinch.hit.col]
   Int_t           Ndata_bb_grinch_hit_pmtnum;
   Double_t        bb_grinch_hit_pmtnum[170];   //[Ndata.bb.grinch.hit.pmtnum]
   Int_t           Ndata_bb_grinch_hit_row;
   Double_t        bb_grinch_hit_row[170];   //[Ndata.bb.grinch.hit.row]
   Int_t           Ndata_bb_grinch_hit_tdc_f;
   Double_t        bb_grinch_hit_tdc_f[170];   //[Ndata.bb.grinch.hit.tdc_f]
   Int_t           Ndata_bb_grinch_hit_tdc_r;
   Double_t        bb_grinch_hit_tdc_r[170];   //[Ndata.bb.grinch.hit.tdc_r]
   Int_t           Ndata_bb_grinch_hit_xhit;
   Double_t        bb_grinch_hit_xhit[170];   //[Ndata.bb.grinch.hit.xhit]
   Int_t           Ndata_bb_grinch_hit_yhit;
   Double_t        bb_grinch_hit_yhit[170];   //[Ndata.bb.grinch.hit.yhit]
   Int_t           Ndata_bb_hodo_bar_tdc_L_le;
   Double_t        bb_hodo_bar_tdc_L_le[27];   //[Ndata.bb.hodo.bar.tdc.L.le]
   Int_t           Ndata_bb_hodo_bar_tdc_L_leW;
   Double_t        bb_hodo_bar_tdc_L_leW[27];   //[Ndata.bb.hodo.bar.tdc.L.leW]
   Int_t           Ndata_bb_hodo_bar_tdc_L_te;
   Double_t        bb_hodo_bar_tdc_L_te[27];   //[Ndata.bb.hodo.bar.tdc.L.te]
   Int_t           Ndata_bb_hodo_bar_tdc_L_teW;
   Double_t        bb_hodo_bar_tdc_L_teW[27];   //[Ndata.bb.hodo.bar.tdc.L.teW]
   Int_t           Ndata_bb_hodo_bar_tdc_L_tot;
   Double_t        bb_hodo_bar_tdc_L_tot[27];   //[Ndata.bb.hodo.bar.tdc.L.tot]
   Int_t           Ndata_bb_hodo_bar_tdc_L_totW;
   Double_t        bb_hodo_bar_tdc_L_totW[27];   //[Ndata.bb.hodo.bar.tdc.L.totW]
   Int_t           Ndata_bb_hodo_bar_tdc_R_le;
   Double_t        bb_hodo_bar_tdc_R_le[27];   //[Ndata.bb.hodo.bar.tdc.R.le]
   Int_t           Ndata_bb_hodo_bar_tdc_R_leW;
   Double_t        bb_hodo_bar_tdc_R_leW[27];   //[Ndata.bb.hodo.bar.tdc.R.leW]
   Int_t           Ndata_bb_hodo_bar_tdc_R_te;
   Double_t        bb_hodo_bar_tdc_R_te[27];   //[Ndata.bb.hodo.bar.tdc.R.te]
   Int_t           Ndata_bb_hodo_bar_tdc_R_teW;
   Double_t        bb_hodo_bar_tdc_R_teW[27];   //[Ndata.bb.hodo.bar.tdc.R.teW]
   Int_t           Ndata_bb_hodo_bar_tdc_R_tot;
   Double_t        bb_hodo_bar_tdc_R_tot[27];   //[Ndata.bb.hodo.bar.tdc.R.tot]
   Int_t           Ndata_bb_hodo_bar_tdc_R_totW;
   Double_t        bb_hodo_bar_tdc_R_totW[27];   //[Ndata.bb.hodo.bar.tdc.R.totW]
   Int_t           Ndata_bb_hodo_bar_tdc_id;
   Double_t        bb_hodo_bar_tdc_id[27];   //[Ndata.bb.hodo.bar.tdc.id]
   Int_t           Ndata_bb_hodo_bar_tdc_meantime;
   Double_t        bb_hodo_bar_tdc_meantime[27];   //[Ndata.bb.hodo.bar.tdc.meantime]
   Int_t           Ndata_bb_hodo_bar_tdc_timediff;
   Double_t        bb_hodo_bar_tdc_timediff[27];   //[Ndata.bb.hodo.bar.tdc.timediff]
   Int_t           Ndata_bb_hodo_bar_tdc_timehitpos;
   Double_t        bb_hodo_bar_tdc_timehitpos[27];   //[Ndata.bb.hodo.bar.tdc.timehitpos]
   Int_t           Ndata_bb_hodo_bar_tdc_vpos;
   Double_t        bb_hodo_bar_tdc_vpos[27];   //[Ndata.bb.hodo.bar.tdc.vpos]
   Int_t           Ndata_bb_hodo_clus_size;
   Double_t        bb_hodo_clus_size[17];   //[Ndata.bb.hodo.clus.size]
   Int_t           Ndata_bb_hodo_clus_tmean;
   Double_t        bb_hodo_clus_tmean[17];   //[Ndata.bb.hodo.clus.tmean]
   Int_t           Ndata_bb_hodo_clus_totmean;
   Double_t        bb_hodo_clus_totmean[17];   //[Ndata.bb.hodo.clus.totmean]
   Int_t           Ndata_bb_hodo_clus_xmean;
   Double_t        bb_hodo_clus_xmean[17];   //[Ndata.bb.hodo.clus.xmean]
   Int_t           Ndata_bb_hodo_clus_ymean;
   Double_t        bb_hodo_clus_ymean[17];   //[Ndata.bb.hodo.clus.ymean]
   Int_t           Ndata_bb_hodo_tdc;
   Double_t        bb_hodo_tdc[180];   //[Ndata.bb.hodo.tdc]
   Int_t           Ndata_bb_hodo_tdc_mult;
   Double_t        bb_hodo_tdc_mult[180];   //[Ndata.bb.hodo.tdc_mult]
   Int_t           Ndata_bb_hodo_tdc_te;
   Double_t        bb_hodo_tdc_te[180];   //[Ndata.bb.hodo.tdc_te]
   Int_t           Ndata_bb_hodo_tdc_tot;
   Double_t        bb_hodo_tdc_tot[180];   //[Ndata.bb.hodo.tdc_tot]
   Int_t           Ndata_bb_hodo_tdccol;
   Double_t        bb_hodo_tdccol[180];   //[Ndata.bb.hodo.tdccol]
   Int_t           Ndata_bb_hodo_tdcelemID;
   Double_t        bb_hodo_tdcelemID[180];   //[Ndata.bb.hodo.tdcelemID]
   Int_t           Ndata_bb_hodo_tdclayer;
   Double_t        bb_hodo_tdclayer[180];   //[Ndata.bb.hodo.tdclayer]
   Int_t           Ndata_bb_hodo_tdcrow;
   Double_t        bb_hodo_tdcrow[180];   //[Ndata.bb.hodo.tdcrow]
   Int_t           Ndata_bb_prob_e;
   Double_t        bb_prob_e[1];   //[Ndata.bb.prob_e]
   Int_t           Ndata_bb_prob_pi;
   Double_t        bb_prob_pi[1];   //[Ndata.bb.prob_pi]
   Int_t           Ndata_bb_ps_a;
   Double_t        bb_ps_a[23];   //[Ndata.bb.ps.a]
   Int_t           Ndata_bb_ps_a_amp;
   Double_t        bb_ps_a_amp[23];   //[Ndata.bb.ps.a_amp]
   Int_t           Ndata_bb_ps_a_amp_p;
   Double_t        bb_ps_a_amp_p[23];   //[Ndata.bb.ps.a_amp_p]
   Int_t           Ndata_bb_ps_a_c;
   Double_t        bb_ps_a_c[23];   //[Ndata.bb.ps.a_c]
   Int_t           Ndata_bb_ps_a_mult;
   Double_t        bb_ps_a_mult[1];   //[Ndata.bb.ps.a_mult]
   Int_t           Ndata_bb_ps_a_p;
   Double_t        bb_ps_a_p[23];   //[Ndata.bb.ps.a_p]
   Int_t           Ndata_bb_ps_a_time;
   Double_t        bb_ps_a_time[23];   //[Ndata.bb.ps.a_time]
   Int_t           Ndata_bb_ps_adccol;
   Double_t        bb_ps_adccol[23];   //[Ndata.bb.ps.adccol]
   Int_t           Ndata_bb_ps_adcelemID;
   Double_t        bb_ps_adcelemID[23];   //[Ndata.bb.ps.adcelemID]
   Int_t           Ndata_bb_ps_adclayer;
   Double_t        bb_ps_adclayer[23];   //[Ndata.bb.ps.adclayer]
   Int_t           Ndata_bb_ps_adcrow;
   Double_t        bb_ps_adcrow[23];   //[Ndata.bb.ps.adcrow]
   Int_t           Ndata_bb_ps_clus_col;
   Double_t        bb_ps_clus_col[7];   //[Ndata.bb.ps.clus.col]
   Int_t           Ndata_bb_ps_clus_e;
   Double_t        bb_ps_clus_e[7];   //[Ndata.bb.ps.clus.e]
   Int_t           Ndata_bb_ps_clus_e_c;
   Double_t        bb_ps_clus_e_c[7];   //[Ndata.bb.ps.clus.e_c]
   Int_t           Ndata_bb_ps_clus_eblk;
   Double_t        bb_ps_clus_eblk[7];   //[Ndata.bb.ps.clus.eblk]
   Int_t           Ndata_bb_ps_clus_eblk_c;
   Double_t        bb_ps_clus_eblk_c[7];   //[Ndata.bb.ps.clus.eblk_c]
   Int_t           Ndata_bb_ps_clus_id;
   Double_t        bb_ps_clus_id[7];   //[Ndata.bb.ps.clus.id]
   Int_t           Ndata_bb_ps_clus_nblk;
   Double_t        bb_ps_clus_nblk[7];   //[Ndata.bb.ps.clus.nblk]
   Int_t           Ndata_bb_ps_clus_row;
   Double_t        bb_ps_clus_row[7];   //[Ndata.bb.ps.clus.row]
   Int_t           Ndata_bb_ps_clus_x;
   Double_t        bb_ps_clus_x[7];   //[Ndata.bb.ps.clus.x]
   Int_t           Ndata_bb_ps_clus_y;
   Double_t        bb_ps_clus_y[7];   //[Ndata.bb.ps.clus.y]
   Int_t           Ndata_bb_ps_clus_blk_col;
   Double_t        bb_ps_clus_blk_col[11];   //[Ndata.bb.ps.clus_blk.col]
   Int_t           Ndata_bb_ps_clus_blk_e;
   Double_t        bb_ps_clus_blk_e[22];   //[Ndata.bb.ps.clus_blk.e]
   Int_t           Ndata_bb_ps_clus_blk_e_c;
   Double_t        bb_ps_clus_blk_e_c[1];   //[Ndata.bb.ps.clus_blk.e_c]
   Int_t           Ndata_bb_ps_clus_blk_id;
   Double_t        bb_ps_clus_blk_id[11];   //[Ndata.bb.ps.clus_blk.id]
   Int_t           Ndata_bb_ps_clus_blk_row;
   Double_t        bb_ps_clus_blk_row[11];   //[Ndata.bb.ps.clus_blk.row]
   Int_t           Ndata_bb_ps_clus_blk_x;
   Double_t        bb_ps_clus_blk_x[11];   //[Ndata.bb.ps.clus_blk.x]
   Int_t           Ndata_bb_ps_clus_blk_y;
   Double_t        bb_ps_clus_blk_y[11];   //[Ndata.bb.ps.clus_blk.y]
   Int_t           Ndata_bb_ps_e_res;
   Double_t        bb_ps_e_res[7];   //[Ndata.bb.ps.e_res]
   Int_t           Ndata_bb_ps_goodblock_atime;
   Double_t        bb_ps_goodblock_atime[12];   //[Ndata.bb.ps.goodblock.atime]
   Int_t           Ndata_bb_ps_goodblock_col;
   Double_t        bb_ps_goodblock_col[12];   //[Ndata.bb.ps.goodblock.col]
   Int_t           Ndata_bb_ps_goodblock_e;
   Double_t        bb_ps_goodblock_e[12];   //[Ndata.bb.ps.goodblock.e]
   Int_t           Ndata_bb_ps_goodblock_id;
   Double_t        bb_ps_goodblock_id[12];   //[Ndata.bb.ps.goodblock.id]
   Int_t           Ndata_bb_ps_goodblock_row;
   Double_t        bb_ps_goodblock_row[12];   //[Ndata.bb.ps.goodblock.row]
   Int_t           Ndata_bb_ps_goodblock_x;
   Double_t        bb_ps_goodblock_x[12];   //[Ndata.bb.ps.goodblock.x]
   Int_t           Ndata_bb_ps_goodblock_y;
   Double_t        bb_ps_goodblock_y[12];   //[Ndata.bb.ps.goodblock.y]
   Int_t           Ndata_bb_ps_ped;
   Double_t        bb_ps_ped[23];   //[Ndata.bb.ps.ped]
   Int_t           Ndata_bb_ps_x_res;
   Double_t        bb_ps_x_res[7];   //[Ndata.bb.ps.x_res]
   Int_t           Ndata_bb_ps_y_res;
   Double_t        bb_ps_y_res[7];   //[Ndata.bb.ps.y_res]
   Int_t           Ndata_bb_sh_a;
   Double_t        bb_sh_a[33];   //[Ndata.bb.sh.a]
   Int_t           Ndata_bb_sh_a_amp;
   Double_t        bb_sh_a_amp[33];   //[Ndata.bb.sh.a_amp]
   Int_t           Ndata_bb_sh_a_amp_p;
   Double_t        bb_sh_a_amp_p[33];   //[Ndata.bb.sh.a_amp_p]
   Int_t           Ndata_bb_sh_a_c;
   Double_t        bb_sh_a_c[33];   //[Ndata.bb.sh.a_c]
   Int_t           Ndata_bb_sh_a_mult;
   Double_t        bb_sh_a_mult[1];   //[Ndata.bb.sh.a_mult]
   Int_t           Ndata_bb_sh_a_p;
   Double_t        bb_sh_a_p[33];   //[Ndata.bb.sh.a_p]
   Int_t           Ndata_bb_sh_a_time;
   Double_t        bb_sh_a_time[33];   //[Ndata.bb.sh.a_time]
   Int_t           Ndata_bb_sh_adccol;
   Double_t        bb_sh_adccol[33];   //[Ndata.bb.sh.adccol]
   Int_t           Ndata_bb_sh_adcelemID;
   Double_t        bb_sh_adcelemID[33];   //[Ndata.bb.sh.adcelemID]
   Int_t           Ndata_bb_sh_adclayer;
   Double_t        bb_sh_adclayer[33];   //[Ndata.bb.sh.adclayer]
   Int_t           Ndata_bb_sh_adcrow;
   Double_t        bb_sh_adcrow[33];   //[Ndata.bb.sh.adcrow]
   Int_t           Ndata_bb_sh_clus_col;
   Double_t        bb_sh_clus_col[7];   //[Ndata.bb.sh.clus.col]
   Int_t           Ndata_bb_sh_clus_e;
   Double_t        bb_sh_clus_e[7];   //[Ndata.bb.sh.clus.e]
   Int_t           Ndata_bb_sh_clus_e_c;
   Double_t        bb_sh_clus_e_c[7];   //[Ndata.bb.sh.clus.e_c]
   Int_t           Ndata_bb_sh_clus_eblk;
   Double_t        bb_sh_clus_eblk[7];   //[Ndata.bb.sh.clus.eblk]
   Int_t           Ndata_bb_sh_clus_eblk_c;
   Double_t        bb_sh_clus_eblk_c[7];   //[Ndata.bb.sh.clus.eblk_c]
   Int_t           Ndata_bb_sh_clus_id;
   Double_t        bb_sh_clus_id[7];   //[Ndata.bb.sh.clus.id]
   Int_t           Ndata_bb_sh_clus_nblk;
   Double_t        bb_sh_clus_nblk[7];   //[Ndata.bb.sh.clus.nblk]
   Int_t           Ndata_bb_sh_clus_row;
   Double_t        bb_sh_clus_row[7];   //[Ndata.bb.sh.clus.row]
   Int_t           Ndata_bb_sh_clus_x;
   Double_t        bb_sh_clus_x[7];   //[Ndata.bb.sh.clus.x]
   Int_t           Ndata_bb_sh_clus_y;
   Double_t        bb_sh_clus_y[7];   //[Ndata.bb.sh.clus.y]
   Int_t           Ndata_bb_sh_clus_blk_col;
   Double_t        bb_sh_clus_blk_col[19];   //[Ndata.bb.sh.clus_blk.col]
   Int_t           Ndata_bb_sh_clus_blk_e;
   Double_t        bb_sh_clus_blk_e[38];   //[Ndata.bb.sh.clus_blk.e]
   Int_t           Ndata_bb_sh_clus_blk_e_c;
   Double_t        bb_sh_clus_blk_e_c[1];   //[Ndata.bb.sh.clus_blk.e_c]
   Int_t           Ndata_bb_sh_clus_blk_id;
   Double_t        bb_sh_clus_blk_id[19];   //[Ndata.bb.sh.clus_blk.id]
   Int_t           Ndata_bb_sh_clus_blk_row;
   Double_t        bb_sh_clus_blk_row[19];   //[Ndata.bb.sh.clus_blk.row]
   Int_t           Ndata_bb_sh_clus_blk_x;
   Double_t        bb_sh_clus_blk_x[19];   //[Ndata.bb.sh.clus_blk.x]
   Int_t           Ndata_bb_sh_clus_blk_y;
   Double_t        bb_sh_clus_blk_y[19];   //[Ndata.bb.sh.clus_blk.y]
   Int_t           Ndata_bb_sh_e_res;
   Double_t        bb_sh_e_res[7];   //[Ndata.bb.sh.e_res]
   Int_t           Ndata_bb_sh_goodblock_atime;
   Double_t        bb_sh_goodblock_atime[22];   //[Ndata.bb.sh.goodblock.atime]
   Int_t           Ndata_bb_sh_goodblock_col;
   Double_t        bb_sh_goodblock_col[22];   //[Ndata.bb.sh.goodblock.col]
   Int_t           Ndata_bb_sh_goodblock_e;
   Double_t        bb_sh_goodblock_e[22];   //[Ndata.bb.sh.goodblock.e]
   Int_t           Ndata_bb_sh_goodblock_id;
   Double_t        bb_sh_goodblock_id[22];   //[Ndata.bb.sh.goodblock.id]
   Int_t           Ndata_bb_sh_goodblock_row;
   Double_t        bb_sh_goodblock_row[22];   //[Ndata.bb.sh.goodblock.row]
   Int_t           Ndata_bb_sh_goodblock_x;
   Double_t        bb_sh_goodblock_x[22];   //[Ndata.bb.sh.goodblock.x]
   Int_t           Ndata_bb_sh_goodblock_y;
   Double_t        bb_sh_goodblock_y[22];   //[Ndata.bb.sh.goodblock.y]
   Int_t           Ndata_bb_sh_ped;
   Double_t        bb_sh_ped[33];   //[Ndata.bb.sh.ped]
   Int_t           Ndata_bb_sh_x_res;
   Double_t        bb_sh_x_res[7];   //[Ndata.bb.sh.x_res]
   Int_t           Ndata_bb_sh_y_res;
   Double_t        bb_sh_y_res[7];   //[Ndata.bb.sh.y_res]
   Int_t           Ndata_bb_tr_beta;
   Double_t        bb_tr_beta[3];   //[Ndata.bb.tr.beta]
   Int_t           Ndata_bb_tr_chi2;
   Double_t        bb_tr_chi2[3];   //[Ndata.bb.tr.chi2]
   Int_t           Ndata_bb_tr_d_ph;
   Double_t        bb_tr_d_ph[3];   //[Ndata.bb.tr.d_ph]
   Int_t           Ndata_bb_tr_d_th;
   Double_t        bb_tr_d_th[3];   //[Ndata.bb.tr.d_th]
   Int_t           Ndata_bb_tr_d_x;
   Double_t        bb_tr_d_x[3];   //[Ndata.bb.tr.d_x]
   Int_t           Ndata_bb_tr_d_y;
   Double_t        bb_tr_d_y[3];   //[Ndata.bb.tr.d_y]
   Int_t           Ndata_bb_tr_dbeta;
   Double_t        bb_tr_dbeta[3];   //[Ndata.bb.tr.dbeta]
   Int_t           Ndata_bb_tr_dtime;
   Double_t        bb_tr_dtime[3];   //[Ndata.bb.tr.dtime]
   Int_t           Ndata_bb_tr_flag;
   Double_t        bb_tr_flag[3];   //[Ndata.bb.tr.flag]
   Int_t           Ndata_bb_tr_ndof;
   Double_t        bb_tr_ndof[3];   //[Ndata.bb.tr.ndof]
   Int_t           Ndata_bb_tr_p;
   Double_t        bb_tr_p[3];   //[Ndata.bb.tr.p]
   Int_t           Ndata_bb_tr_pathl;
   Double_t        bb_tr_pathl[3];   //[Ndata.bb.tr.pathl]
   Int_t           Ndata_bb_tr_ph;
   Double_t        bb_tr_ph[3];   //[Ndata.bb.tr.ph]
   Int_t           Ndata_bb_tr_px;
   Double_t        bb_tr_px[3];   //[Ndata.bb.tr.px]
   Int_t           Ndata_bb_tr_py;
   Double_t        bb_tr_py[3];   //[Ndata.bb.tr.py]
   Int_t           Ndata_bb_tr_pz;
   Double_t        bb_tr_pz[3];   //[Ndata.bb.tr.pz]
   Int_t           Ndata_bb_tr_r_ph;
   Double_t        bb_tr_r_ph[3];   //[Ndata.bb.tr.r_ph]
   Int_t           Ndata_bb_tr_r_th;
   Double_t        bb_tr_r_th[3];   //[Ndata.bb.tr.r_th]
   Int_t           Ndata_bb_tr_r_x;
   Double_t        bb_tr_r_x[3];   //[Ndata.bb.tr.r_x]
   Int_t           Ndata_bb_tr_r_y;
   Double_t        bb_tr_r_y[3];   //[Ndata.bb.tr.r_y]
   Int_t           Ndata_bb_tr_tg_dp;
   Double_t        bb_tr_tg_dp[3];   //[Ndata.bb.tr.tg_dp]
   Int_t           Ndata_bb_tr_tg_ph;
   Double_t        bb_tr_tg_ph[3];   //[Ndata.bb.tr.tg_ph]
   Int_t           Ndata_bb_tr_tg_th;
   Double_t        bb_tr_tg_th[3];   //[Ndata.bb.tr.tg_th]
   Int_t           Ndata_bb_tr_tg_y;
   Double_t        bb_tr_tg_y[3];   //[Ndata.bb.tr.tg_y]
   Int_t           Ndata_bb_tr_th;
   Double_t        bb_tr_th[3];   //[Ndata.bb.tr.th]
   Int_t           Ndata_bb_tr_time;
   Double_t        bb_tr_time[3];   //[Ndata.bb.tr.time]
   Int_t           Ndata_bb_tr_vx;
   Double_t        bb_tr_vx[3];   //[Ndata.bb.tr.vx]
   Int_t           Ndata_bb_tr_vy;
   Double_t        bb_tr_vy[3];   //[Ndata.bb.tr.vy]
   Int_t           Ndata_bb_tr_vz;
   Double_t        bb_tr_vz[3];   //[Ndata.bb.tr.vz]
   Int_t           Ndata_bb_tr_x;
   Double_t        bb_tr_x[3];   //[Ndata.bb.tr.x]
   Int_t           Ndata_bb_tr_y;
   Double_t        bb_tr_y[3];   //[Ndata.bb.tr.y]
   Int_t           Ndata_bb_x_bcp;
   Double_t        bb_x_bcp[1];   //[Ndata.bb.x_bcp]
   Int_t           Ndata_bb_x_fcp;
   Double_t        bb_x_fcp[1];   //[Ndata.bb.x_fcp]
   Int_t           Ndata_bb_y_bcp;
   Double_t        bb_y_bcp[1];   //[Ndata.bb.y_bcp]
   Int_t           Ndata_bb_y_fcp;
   Double_t        bb_y_fcp[1];   //[Ndata.bb.y_fcp]
   Int_t           Ndata_sbs_hcal_a;
   Double_t        sbs_hcal_a[67];   //[Ndata.sbs.hcal.a]
   Int_t           Ndata_sbs_hcal_a_amp;
   Double_t        sbs_hcal_a_amp[67];   //[Ndata.sbs.hcal.a_amp]
   Int_t           Ndata_sbs_hcal_a_amp_p;
   Double_t        sbs_hcal_a_amp_p[67];   //[Ndata.sbs.hcal.a_amp_p]
   Int_t           Ndata_sbs_hcal_a_c;
   Double_t        sbs_hcal_a_c[67];   //[Ndata.sbs.hcal.a_c]
   Int_t           Ndata_sbs_hcal_a_mult;
   Double_t        sbs_hcal_a_mult[1];   //[Ndata.sbs.hcal.a_mult]
   Int_t           Ndata_sbs_hcal_a_p;
   Double_t        sbs_hcal_a_p[67];   //[Ndata.sbs.hcal.a_p]
   Int_t           Ndata_sbs_hcal_a_time;
   Double_t        sbs_hcal_a_time[67];   //[Ndata.sbs.hcal.a_time]
   Int_t           Ndata_sbs_hcal_adccol;
   Double_t        sbs_hcal_adccol[67];   //[Ndata.sbs.hcal.adccol]
   Int_t           Ndata_sbs_hcal_adcelemID;
   Double_t        sbs_hcal_adcelemID[67];   //[Ndata.sbs.hcal.adcelemID]
   Int_t           Ndata_sbs_hcal_adclayer;
   Double_t        sbs_hcal_adclayer[67];   //[Ndata.sbs.hcal.adclayer]
   Int_t           Ndata_sbs_hcal_adcrow;
   Double_t        sbs_hcal_adcrow[67];   //[Ndata.sbs.hcal.adcrow]
   Int_t           Ndata_sbs_hcal_clus_col;
   Double_t        sbs_hcal_clus_col[1];   //[Ndata.sbs.hcal.clus.col]
   Int_t           Ndata_sbs_hcal_clus_e;
   Double_t        sbs_hcal_clus_e[1];   //[Ndata.sbs.hcal.clus.e]
   Int_t           Ndata_sbs_hcal_clus_e_c;
   Double_t        sbs_hcal_clus_e_c[1];   //[Ndata.sbs.hcal.clus.e_c]
   Int_t           Ndata_sbs_hcal_clus_eblk;
   Double_t        sbs_hcal_clus_eblk[1];   //[Ndata.sbs.hcal.clus.eblk]
   Int_t           Ndata_sbs_hcal_clus_eblk_c;
   Double_t        sbs_hcal_clus_eblk_c[1];   //[Ndata.sbs.hcal.clus.eblk_c]
   Int_t           Ndata_sbs_hcal_clus_id;
   Double_t        sbs_hcal_clus_id[1];   //[Ndata.sbs.hcal.clus.id]
   Int_t           Ndata_sbs_hcal_clus_nblk;
   Double_t        sbs_hcal_clus_nblk[1];   //[Ndata.sbs.hcal.clus.nblk]
   Int_t           Ndata_sbs_hcal_clus_row;
   Double_t        sbs_hcal_clus_row[1];   //[Ndata.sbs.hcal.clus.row]
   Int_t           Ndata_sbs_hcal_clus_x;
   Double_t        sbs_hcal_clus_x[1];   //[Ndata.sbs.hcal.clus.x]
   Int_t           Ndata_sbs_hcal_clus_y;
   Double_t        sbs_hcal_clus_y[1];   //[Ndata.sbs.hcal.clus.y]
   Int_t           Ndata_sbs_hcal_clus_blk_col;
   Double_t        sbs_hcal_clus_blk_col[25];   //[Ndata.sbs.hcal.clus_blk.col]
   Int_t           Ndata_sbs_hcal_clus_blk_e;
   Double_t        sbs_hcal_clus_blk_e[50];   //[Ndata.sbs.hcal.clus_blk.e]
   Int_t           Ndata_sbs_hcal_clus_blk_e_c;
   Double_t        sbs_hcal_clus_blk_e_c[1];   //[Ndata.sbs.hcal.clus_blk.e_c]
   Int_t           Ndata_sbs_hcal_clus_blk_id;
   Double_t        sbs_hcal_clus_blk_id[25];   //[Ndata.sbs.hcal.clus_blk.id]
   Int_t           Ndata_sbs_hcal_clus_blk_row;
   Double_t        sbs_hcal_clus_blk_row[25];   //[Ndata.sbs.hcal.clus_blk.row]
   Int_t           Ndata_sbs_hcal_clus_blk_x;
   Double_t        sbs_hcal_clus_blk_x[25];   //[Ndata.sbs.hcal.clus_blk.x]
   Int_t           Ndata_sbs_hcal_clus_blk_y;
   Double_t        sbs_hcal_clus_blk_y[25];   //[Ndata.sbs.hcal.clus_blk.y]
   Int_t           Ndata_sbs_hcal_goodblock_atime;
   Double_t        sbs_hcal_goodblock_atime[25];   //[Ndata.sbs.hcal.goodblock.atime]
   Int_t           Ndata_sbs_hcal_goodblock_col;
   Double_t        sbs_hcal_goodblock_col[25];   //[Ndata.sbs.hcal.goodblock.col]
   Int_t           Ndata_sbs_hcal_goodblock_e;
   Double_t        sbs_hcal_goodblock_e[25];   //[Ndata.sbs.hcal.goodblock.e]
   Int_t           Ndata_sbs_hcal_goodblock_id;
   Double_t        sbs_hcal_goodblock_id[25];   //[Ndata.sbs.hcal.goodblock.id]
   Int_t           Ndata_sbs_hcal_goodblock_row;
   Double_t        sbs_hcal_goodblock_row[25];   //[Ndata.sbs.hcal.goodblock.row]
   Int_t           Ndata_sbs_hcal_goodblock_x;
   Double_t        sbs_hcal_goodblock_x[25];   //[Ndata.sbs.hcal.goodblock.x]
   Int_t           Ndata_sbs_hcal_goodblock_y;
   Double_t        sbs_hcal_goodblock_y[25];   //[Ndata.sbs.hcal.goodblock.y]
   Int_t           Ndata_sbs_hcal_ped;
   Double_t        sbs_hcal_ped[67];   //[Ndata.sbs.hcal.ped]
   Int_t           Ndata_sbs_hcal_tdc;
   Double_t        sbs_hcal_tdc[288];   //[Ndata.sbs.hcal.tdc]
   Int_t           Ndata_sbs_hcal_tdc_mult;
   Double_t        sbs_hcal_tdc_mult[288];   //[Ndata.sbs.hcal.tdc_mult]
   Int_t           Ndata_sbs_hcal_tdccol;
   Double_t        sbs_hcal_tdccol[288];   //[Ndata.sbs.hcal.tdccol]
   Int_t           Ndata_sbs_hcal_tdcelemID;
   Double_t        sbs_hcal_tdcelemID[288];   //[Ndata.sbs.hcal.tdcelemID]
   Int_t           Ndata_sbs_hcal_tdclayer;
   Double_t        sbs_hcal_tdclayer[288];   //[Ndata.sbs.hcal.tdclayer]
   Int_t           Ndata_sbs_hcal_tdcrow;
   Double_t        sbs_hcal_tdcrow[288];   //[Ndata.sbs.hcal.tdcrow]
   Int_t           Ndata_sbs_tr_beta;
   Double_t        sbs_tr_beta[1];   //[Ndata.sbs.tr.beta]
   Int_t           Ndata_sbs_tr_chi2;
   Double_t        sbs_tr_chi2[1];   //[Ndata.sbs.tr.chi2]
   Int_t           Ndata_sbs_tr_d_ph;
   Double_t        sbs_tr_d_ph[1];   //[Ndata.sbs.tr.d_ph]
   Int_t           Ndata_sbs_tr_d_th;
   Double_t        sbs_tr_d_th[1];   //[Ndata.sbs.tr.d_th]
   Int_t           Ndata_sbs_tr_d_x;
   Double_t        sbs_tr_d_x[1];   //[Ndata.sbs.tr.d_x]
   Int_t           Ndata_sbs_tr_d_y;
   Double_t        sbs_tr_d_y[1];   //[Ndata.sbs.tr.d_y]
   Int_t           Ndata_sbs_tr_dbeta;
   Double_t        sbs_tr_dbeta[1];   //[Ndata.sbs.tr.dbeta]
   Int_t           Ndata_sbs_tr_dtime;
   Double_t        sbs_tr_dtime[1];   //[Ndata.sbs.tr.dtime]
   Int_t           Ndata_sbs_tr_flag;
   Double_t        sbs_tr_flag[1];   //[Ndata.sbs.tr.flag]
   Int_t           Ndata_sbs_tr_ndof;
   Double_t        sbs_tr_ndof[1];   //[Ndata.sbs.tr.ndof]
   Int_t           Ndata_sbs_tr_p;
   Double_t        sbs_tr_p[1];   //[Ndata.sbs.tr.p]
   Int_t           Ndata_sbs_tr_pathl;
   Double_t        sbs_tr_pathl[1];   //[Ndata.sbs.tr.pathl]
   Int_t           Ndata_sbs_tr_ph;
   Double_t        sbs_tr_ph[1];   //[Ndata.sbs.tr.ph]
   Int_t           Ndata_sbs_tr_px;
   Double_t        sbs_tr_px[1];   //[Ndata.sbs.tr.px]
   Int_t           Ndata_sbs_tr_py;
   Double_t        sbs_tr_py[1];   //[Ndata.sbs.tr.py]
   Int_t           Ndata_sbs_tr_pz;
   Double_t        sbs_tr_pz[1];   //[Ndata.sbs.tr.pz]
   Int_t           Ndata_sbs_tr_r_ph;
   Double_t        sbs_tr_r_ph[1];   //[Ndata.sbs.tr.r_ph]
   Int_t           Ndata_sbs_tr_r_th;
   Double_t        sbs_tr_r_th[1];   //[Ndata.sbs.tr.r_th]
   Int_t           Ndata_sbs_tr_r_x;
   Double_t        sbs_tr_r_x[1];   //[Ndata.sbs.tr.r_x]
   Int_t           Ndata_sbs_tr_r_y;
   Double_t        sbs_tr_r_y[1];   //[Ndata.sbs.tr.r_y]
   Int_t           Ndata_sbs_tr_tg_dp;
   Double_t        sbs_tr_tg_dp[1];   //[Ndata.sbs.tr.tg_dp]
   Int_t           Ndata_sbs_tr_tg_ph;
   Double_t        sbs_tr_tg_ph[1];   //[Ndata.sbs.tr.tg_ph]
   Int_t           Ndata_sbs_tr_tg_th;
   Double_t        sbs_tr_tg_th[1];   //[Ndata.sbs.tr.tg_th]
   Int_t           Ndata_sbs_tr_tg_y;
   Double_t        sbs_tr_tg_y[1];   //[Ndata.sbs.tr.tg_y]
   Int_t           Ndata_sbs_tr_th;
   Double_t        sbs_tr_th[1];   //[Ndata.sbs.tr.th]
   Int_t           Ndata_sbs_tr_time;
   Double_t        sbs_tr_time[1];   //[Ndata.sbs.tr.time]
   Int_t           Ndata_sbs_tr_vx;
   Double_t        sbs_tr_vx[1];   //[Ndata.sbs.tr.vx]
   Int_t           Ndata_sbs_tr_vy;
   Double_t        sbs_tr_vy[1];   //[Ndata.sbs.tr.vy]
   Int_t           Ndata_sbs_tr_vz;
   Double_t        sbs_tr_vz[1];   //[Ndata.sbs.tr.vz]
   Int_t           Ndata_sbs_tr_x;
   Double_t        sbs_tr_x[1];   //[Ndata.sbs.tr.x]
   Int_t           Ndata_sbs_tr_y;
   Double_t        sbs_tr_y[1];   //[Ndata.sbs.tr.y]
   Double_t        MC_hit_n;
   Double_t        MC_mc_px;
   Double_t        MC_mc_py;
   Double_t        MC_mc_pz;
   Double_t        MC_mc_vx;
   Double_t        MC_mc_vy;
   Double_t        MC_mc_vz;
   Double_t        MC_nbbgemhits;
   Double_t        MC_nbbtracks;
   Double_t        MC_pt_n;
   Double_t        MC_tr_n;
   Double_t        MC_weight;
   Double_t        bb_gem_hit_ngoodhits;
   Double_t        bb_gem_nlayershit;
   Double_t        bb_gem_nlayershitu;
   Double_t        bb_gem_nlayershituv;
   Double_t        bb_gem_nlayershitv;
   Double_t        bb_gem_track_besttrack;
   Double_t        bb_gem_track_ntrack;
   Double_t        bb_grinch_nclus;
   Double_t        bb_grinch_nhits;
   Double_t        bb_hodo_bar_ngoodbars;
   Double_t        bb_hodo_nclus;
   Double_t        bb_hodo_ngoodADChits;
   Double_t        bb_hodo_ngoodTDChits;
   Double_t        bb_hodo_nhits;
   Double_t        bb_hodo_nrefhits;
   Double_t        bb_ps_colblk;
   Double_t        bb_ps_e;
   Double_t        bb_ps_e_c;
   Double_t        bb_ps_e_m_res;
   Double_t        bb_ps_eblk;
   Double_t        bb_ps_eblk_c;
   Double_t        bb_ps_idblk;
   Double_t        bb_ps_nblk;
   Double_t        bb_ps_nclus;
   Double_t        bb_ps_ngoodADChits;
   Double_t        bb_ps_ngoodTDChits;
   Double_t        bb_ps_nhits;
   Double_t        bb_ps_nrefhits;
   Double_t        bb_ps_rowblk;
   Double_t        bb_ps_x;
   Double_t        bb_ps_x_m_res;
   Double_t        bb_ps_y;
   Double_t        bb_ps_y_m_res;
   Double_t        bb_sh_colblk;
   Double_t        bb_sh_e;
   Double_t        bb_sh_e_c;
   Double_t        bb_sh_e_m_res;
   Double_t        bb_sh_eblk;
   Double_t        bb_sh_eblk_c;
   Double_t        bb_sh_idblk;
   Double_t        bb_sh_nblk;
   Double_t        bb_sh_nclus;
   Double_t        bb_sh_ngoodADChits;
   Double_t        bb_sh_ngoodTDChits;
   Double_t        bb_sh_nhits;
   Double_t        bb_sh_nrefhits;
   Double_t        bb_sh_rowblk;
   Double_t        bb_sh_x;
   Double_t        bb_sh_x_m_res;
   Double_t        bb_sh_y;
   Double_t        bb_sh_y_m_res;
   Double_t        bb_tr_n;
   Double_t        sbs_hcal_colblk;
   Double_t        sbs_hcal_e;
   Double_t        sbs_hcal_e_c;
   Double_t        sbs_hcal_eblk;
   Double_t        sbs_hcal_eblk_c;
   Double_t        sbs_hcal_idblk;
   Double_t        sbs_hcal_ledbit;
   Double_t        sbs_hcal_ledcount;
   Double_t        sbs_hcal_nblk;
   Double_t        sbs_hcal_nclus;
   Double_t        sbs_hcal_ngoodADChits;
   Double_t        sbs_hcal_ngoodTDChits;
   Double_t        sbs_hcal_nhits;
   Double_t        sbs_hcal_nrefhits;
   Double_t        sbs_hcal_rowblk;
   Double_t        sbs_hcal_x;
   Double_t        sbs_hcal_y;
   Double_t        sbs_status;
   Double_t        sbs_tr_n;
 //THaEvent        *Event_Branch;
   ULong64_t       fEvtHdr_fEvtTime;
   UInt_t          fEvtHdr_fEvtNum;
   UInt_t          fEvtHdr_fEvtType;
   UInt_t          fEvtHdr_fEvtLen;
   Int_t           fEvtHdr_fHelicity;
   Int_t           fEvtHdr_fTargetPol;
   UInt_t          fEvtHdr_fRun;

   // List of branches
   TBranch        *b_Ndata_MC_bbgemhit_edep;   //!
   TBranch        *b_MC_bbgemhit_edep;   //!
   TBranch        *b_Ndata_MC_bbgemhit_mid;   //!
   TBranch        *b_MC_bbgemhit_mid;   //!
   TBranch        *b_Ndata_MC_bbgemhit_pid;   //!
   TBranch        *b_MC_bbgemhit_pid;   //!
   TBranch        *b_Ndata_MC_bbgemhit_plane;   //!
   TBranch        *b_MC_bbgemhit_plane;   //!
   TBranch        *b_Ndata_MC_bbgemhit_tid;   //!
   TBranch        *b_MC_bbgemhit_tid;   //!
   TBranch        *b_Ndata_MC_bbgemhit_x;   //!
   TBranch        *b_MC_bbgemhit_x;   //!
   TBranch        *b_Ndata_MC_bbgemhit_y;   //!
   TBranch        *b_MC_bbgemhit_y;   //!
   TBranch        *b_Ndata_MC_bbtrack_dx;   //!
   TBranch        *b_MC_bbtrack_dx;   //!
   TBranch        *b_Ndata_MC_bbtrack_dy;   //!
   TBranch        *b_MC_bbtrack_dy;   //!
   TBranch        *b_Ndata_MC_bbtrack_mid;   //!
   TBranch        *b_MC_bbtrack_mid;   //!
   TBranch        *b_Ndata_MC_bbtrack_nhits;   //!
   TBranch        *b_MC_bbtrack_nhits;   //!
   TBranch        *b_Ndata_MC_bbtrack_p;   //!
   TBranch        *b_MC_bbtrack_p;   //!
   TBranch        *b_Ndata_MC_bbtrack_pid;   //!
   TBranch        *b_MC_bbtrack_pid;   //!
   TBranch        *b_Ndata_MC_bbtrack_tid;   //!
   TBranch        *b_MC_bbtrack_tid;   //!
   TBranch        *b_Ndata_MC_bbtrack_x;   //!
   TBranch        *b_MC_bbtrack_x;   //!
   TBranch        *b_Ndata_MC_bbtrack_y;   //!
   TBranch        *b_MC_bbtrack_y;   //!
   TBranch        *b_Ndata_MC_pt_clustsz;   //!
   TBranch        *b_MC_pt_clustsz;   //!
   TBranch        *b_Ndata_MC_pt_deflect;   //!
   TBranch        *b_MC_pt_deflect;   //!
   TBranch        *b_Ndata_MC_pt_deltaE;   //!
   TBranch        *b_MC_pt_deltaE;   //!
   TBranch        *b_Ndata_MC_pt_hitres;   //!
   TBranch        *b_MC_pt_hitres;   //!
   TBranch        *b_Ndata_MC_pt_nfound;   //!
   TBranch        *b_MC_pt_nfound;   //!
   TBranch        *b_Ndata_MC_pt_p;   //!
   TBranch        *b_MC_pt_p;   //!
   TBranch        *b_Ndata_MC_pt_ph;   //!
   TBranch        *b_MC_pt_ph;   //!
   TBranch        *b_Ndata_MC_pt_phdir;   //!
   TBranch        *b_MC_pt_phdir;   //!
   TBranch        *b_Ndata_MC_pt_phi;   //!
   TBranch        *b_MC_pt_phi;   //!
   TBranch        *b_Ndata_MC_pt_plane;   //!
   TBranch        *b_MC_pt_plane;   //!
   TBranch        *b_Ndata_MC_pt_r;   //!
   TBranch        *b_MC_pt_r;   //!
   TBranch        *b_Ndata_MC_pt_status;   //!
   TBranch        *b_MC_pt_status;   //!
   TBranch        *b_Ndata_MC_pt_th;   //!
   TBranch        *b_MC_pt_th;   //!
   TBranch        *b_Ndata_MC_pt_thdir;   //!
   TBranch        *b_MC_pt_thdir;   //!
   TBranch        *b_Ndata_MC_pt_theta;   //!
   TBranch        *b_MC_pt_theta;   //!
   TBranch        *b_Ndata_MC_pt_time;   //!
   TBranch        *b_MC_pt_time;   //!
   TBranch        *b_Ndata_MC_pt_tof;   //!
   TBranch        *b_MC_pt_tof;   //!
   TBranch        *b_Ndata_MC_pt_trkres;   //!
   TBranch        *b_MC_pt_trkres;   //!
   TBranch        *b_Ndata_MC_pt_type;   //!
   TBranch        *b_MC_pt_type;   //!
   TBranch        *b_Ndata_MC_pt_x;   //!
   TBranch        *b_MC_pt_x;   //!
   TBranch        *b_Ndata_MC_pt_y;   //!
   TBranch        *b_MC_pt_y;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCU;   //!
   TBranch        *b_bb_gem_hit_ADCU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCV;   //!
   TBranch        *b_bb_gem_hit_ADCV;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCasym;   //!
   TBranch        *b_bb_gem_hit_ADCasym;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCavg;   //!
   TBranch        *b_bb_gem_hit_ADCavg;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampU;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampUclust;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampUclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampV;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampV;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxsampVclust;   //!
   TBranch        *b_bb_gem_hit_ADCmaxsampVclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxstripU;   //!
   TBranch        *b_bb_gem_hit_ADCmaxstripU;   //!
   TBranch        *b_Ndata_bb_gem_hit_ADCmaxstripV;   //!
   TBranch        *b_bb_gem_hit_ADCmaxstripV;   //!
   TBranch        *b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_bb_gem_hit_BUILD_ALL_SAMPLES_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_bb_gem_hit_BUILD_ALL_SAMPLES_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_CM_OR_U;   //!
   TBranch        *b_bb_gem_hit_CM_OR_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_CM_OR_V;   //!
   TBranch        *b_bb_gem_hit_CM_OR_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_ENABLE_CM_U;   //!
   TBranch        *b_bb_gem_hit_ENABLE_CM_U;   //!
   TBranch        *b_Ndata_bb_gem_hit_ENABLE_CM_V;   //!
   TBranch        *b_bb_gem_hit_ENABLE_CM_V;   //!
   TBranch        *b_Ndata_bb_gem_hit_Tavg;   //!
   TBranch        *b_bb_gem_hit_Tavg;   //!
   TBranch        *b_Ndata_bb_gem_hit_Utime;   //!
   TBranch        *b_bb_gem_hit_Utime;   //!
   TBranch        *b_Ndata_bb_gem_hit_UtimeMaxStrip;   //!
   TBranch        *b_bb_gem_hit_UtimeMaxStrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_Vtime;   //!
   TBranch        *b_bb_gem_hit_Vtime;   //!
   TBranch        *b_Ndata_bb_gem_hit_VtimeMaxStrip;   //!
   TBranch        *b_bb_gem_hit_VtimeMaxStrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_ccor_clust;   //!
   TBranch        *b_bb_gem_hit_ccor_clust;   //!
   TBranch        *b_Ndata_bb_gem_hit_ccor_strip;   //!
   TBranch        *b_bb_gem_hit_ccor_strip;   //!
   TBranch        *b_Ndata_bb_gem_hit_deltat;   //!
   TBranch        *b_bb_gem_hit_deltat;   //!
   TBranch        *b_Ndata_bb_gem_hit_eresidu;   //!
   TBranch        *b_bb_gem_hit_eresidu;   //!
   TBranch        *b_Ndata_bb_gem_hit_eresidv;   //!
   TBranch        *b_bb_gem_hit_eresidv;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxUclust;   //!
   TBranch        *b_bb_gem_hit_isampmaxUclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxUstrip;   //!
   TBranch        *b_bb_gem_hit_isampmaxUstrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxVclust;   //!
   TBranch        *b_bb_gem_hit_isampmaxVclust;   //!
   TBranch        *b_Ndata_bb_gem_hit_isampmaxVstrip;   //!
   TBranch        *b_bb_gem_hit_isampmaxVstrip;   //!
   TBranch        *b_Ndata_bb_gem_hit_layer;   //!
   TBranch        *b_bb_gem_hit_layer;   //!
   TBranch        *b_Ndata_bb_gem_hit_module;   //!
   TBranch        *b_bb_gem_hit_module;   //!
   TBranch        *b_Ndata_bb_gem_hit_nstripu;   //!
   TBranch        *b_bb_gem_hit_nstripu;   //!
   TBranch        *b_Ndata_bb_gem_hit_nstripv;   //!
   TBranch        *b_bb_gem_hit_nstripv;   //!
   TBranch        *b_Ndata_bb_gem_hit_residu;   //!
   TBranch        *b_bb_gem_hit_residu;   //!
   TBranch        *b_Ndata_bb_gem_hit_residv;   //!
   TBranch        *b_bb_gem_hit_residv;   //!
   TBranch        *b_Ndata_bb_gem_hit_trackindex;   //!
   TBranch        *b_bb_gem_hit_trackindex;   //!
   TBranch        *b_Ndata_bb_gem_hit_u;   //!
   TBranch        *b_bb_gem_hit_u;   //!
   TBranch        *b_Ndata_bb_gem_hit_umoment;   //!
   TBranch        *b_bb_gem_hit_umoment;   //!
   TBranch        *b_Ndata_bb_gem_hit_usigma;   //!
   TBranch        *b_bb_gem_hit_usigma;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustriphi;   //!
   TBranch        *b_bb_gem_hit_ustriphi;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustriplo;   //!
   TBranch        *b_bb_gem_hit_ustriplo;   //!
   TBranch        *b_Ndata_bb_gem_hit_ustripmax;   //!
   TBranch        *b_bb_gem_hit_ustripmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_v;   //!
   TBranch        *b_bb_gem_hit_v;   //!
   TBranch        *b_Ndata_bb_gem_hit_vmoment;   //!
   TBranch        *b_bb_gem_hit_vmoment;   //!
   TBranch        *b_Ndata_bb_gem_hit_vsigma;   //!
   TBranch        *b_bb_gem_hit_vsigma;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstriphi;   //!
   TBranch        *b_bb_gem_hit_vstriphi;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstriplo;   //!
   TBranch        *b_bb_gem_hit_vstriplo;   //!
   TBranch        *b_Ndata_bb_gem_hit_vstripmax;   //!
   TBranch        *b_bb_gem_hit_vstripmax;   //!
   TBranch        *b_Ndata_bb_gem_hit_xglobal;   //!
   TBranch        *b_bb_gem_hit_xglobal;   //!
   TBranch        *b_Ndata_bb_gem_hit_xlocal;   //!
   TBranch        *b_bb_gem_hit_xlocal;   //!
   TBranch        *b_Ndata_bb_gem_hit_yglobal;   //!
   TBranch        *b_bb_gem_hit_yglobal;   //!
   TBranch        *b_Ndata_bb_gem_hit_ylocal;   //!
   TBranch        *b_bb_gem_hit_ylocal;   //!
   TBranch        *b_Ndata_bb_gem_hit_zglobal;   //!
   TBranch        *b_bb_gem_hit_zglobal;   //!
   TBranch        *b_Ndata_bb_gem_n2Dhit_layer;   //!
   TBranch        *b_bb_gem_n2Dhit_layer;   //!
   TBranch        *b_Ndata_bb_gem_nclustu_layer;   //!
   TBranch        *b_bb_gem_nclustu_layer;   //!
   TBranch        *b_Ndata_bb_gem_nclustv_layer;   //!
   TBranch        *b_bb_gem_nclustv_layer;   //!
   TBranch        *b_Ndata_bb_gem_nstripsu_layer;   //!
   TBranch        *b_bb_gem_nstripsu_layer;   //!
   TBranch        *b_Ndata_bb_gem_nstripsv_layer;   //!
   TBranch        *b_bb_gem_nstripsv_layer;   //!
   TBranch        *b_Ndata_bb_gem_track_chi2ndf;   //!
   TBranch        *b_bb_gem_track_chi2ndf;   //!
   TBranch        *b_Ndata_bb_gem_track_nhits;   //!
   TBranch        *b_bb_gem_track_nhits;   //!
   TBranch        *b_Ndata_bb_gem_track_x;   //!
   TBranch        *b_bb_gem_track_x;   //!
   TBranch        *b_Ndata_bb_gem_track_xp;   //!
   TBranch        *b_bb_gem_track_xp;   //!
   TBranch        *b_Ndata_bb_gem_track_y;   //!
   TBranch        *b_bb_gem_track_y;   //!
   TBranch        *b_Ndata_bb_gem_track_yp;   //!
   TBranch        *b_bb_gem_track_yp;   //!
   TBranch        *b_Ndata_bb_grinch_clus_adc;   //!
   TBranch        *b_bb_grinch_clus_adc;   //!
   TBranch        *b_Ndata_bb_grinch_clus_size;   //!
   TBranch        *b_bb_grinch_clus_size;   //!
   TBranch        *b_Ndata_bb_grinch_clus_tf_mean;   //!
   TBranch        *b_bb_grinch_clus_tf_mean;   //!
   TBranch        *b_Ndata_bb_grinch_clus_tr_mean;   //!
   TBranch        *b_bb_grinch_clus_tr_mean;   //!
   TBranch        *b_Ndata_bb_grinch_clus_x_mean;   //!
   TBranch        *b_bb_grinch_clus_x_mean;   //!
   TBranch        *b_Ndata_bb_grinch_clus_y_mean;   //!
   TBranch        *b_bb_grinch_clus_y_mean;   //!
   TBranch        *b_Ndata_bb_grinch_hit_adc;   //!
   TBranch        *b_bb_grinch_hit_adc;   //!
   TBranch        *b_Ndata_bb_grinch_hit_col;   //!
   TBranch        *b_bb_grinch_hit_col;   //!
   TBranch        *b_Ndata_bb_grinch_hit_pmtnum;   //!
   TBranch        *b_bb_grinch_hit_pmtnum;   //!
   TBranch        *b_Ndata_bb_grinch_hit_row;   //!
   TBranch        *b_bb_grinch_hit_row;   //!
   TBranch        *b_Ndata_bb_grinch_hit_tdc_f;   //!
   TBranch        *b_bb_grinch_hit_tdc_f;   //!
   TBranch        *b_Ndata_bb_grinch_hit_tdc_r;   //!
   TBranch        *b_bb_grinch_hit_tdc_r;   //!
   TBranch        *b_Ndata_bb_grinch_hit_xhit;   //!
   TBranch        *b_bb_grinch_hit_xhit;   //!
   TBranch        *b_Ndata_bb_grinch_hit_yhit;   //!
   TBranch        *b_bb_grinch_hit_yhit;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_L_le;   //!
   TBranch        *b_bb_hodo_bar_tdc_L_le;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_L_leW;   //!
   TBranch        *b_bb_hodo_bar_tdc_L_leW;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_L_te;   //!
   TBranch        *b_bb_hodo_bar_tdc_L_te;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_L_teW;   //!
   TBranch        *b_bb_hodo_bar_tdc_L_teW;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_L_tot;   //!
   TBranch        *b_bb_hodo_bar_tdc_L_tot;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_L_totW;   //!
   TBranch        *b_bb_hodo_bar_tdc_L_totW;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_R_le;   //!
   TBranch        *b_bb_hodo_bar_tdc_R_le;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_R_leW;   //!
   TBranch        *b_bb_hodo_bar_tdc_R_leW;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_R_te;   //!
   TBranch        *b_bb_hodo_bar_tdc_R_te;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_R_teW;   //!
   TBranch        *b_bb_hodo_bar_tdc_R_teW;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_R_tot;   //!
   TBranch        *b_bb_hodo_bar_tdc_R_tot;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_R_totW;   //!
   TBranch        *b_bb_hodo_bar_tdc_R_totW;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_id;   //!
   TBranch        *b_bb_hodo_bar_tdc_id;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_meantime;   //!
   TBranch        *b_bb_hodo_bar_tdc_meantime;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_timediff;   //!
   TBranch        *b_bb_hodo_bar_tdc_timediff;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_timehitpos;   //!
   TBranch        *b_bb_hodo_bar_tdc_timehitpos;   //!
   TBranch        *b_Ndata_bb_hodo_bar_tdc_vpos;   //!
   TBranch        *b_bb_hodo_bar_tdc_vpos;   //!
   TBranch        *b_Ndata_bb_hodo_clus_size;   //!
   TBranch        *b_bb_hodo_clus_size;   //!
   TBranch        *b_Ndata_bb_hodo_clus_tmean;   //!
   TBranch        *b_bb_hodo_clus_tmean;   //!
   TBranch        *b_Ndata_bb_hodo_clus_totmean;   //!
   TBranch        *b_bb_hodo_clus_totmean;   //!
   TBranch        *b_Ndata_bb_hodo_clus_xmean;   //!
   TBranch        *b_bb_hodo_clus_xmean;   //!
   TBranch        *b_Ndata_bb_hodo_clus_ymean;   //!
   TBranch        *b_bb_hodo_clus_ymean;   //!
   TBranch        *b_Ndata_bb_hodo_tdc;   //!
   TBranch        *b_bb_hodo_tdc;   //!
   TBranch        *b_Ndata_bb_hodo_tdc_mult;   //!
   TBranch        *b_bb_hodo_tdc_mult;   //!
   TBranch        *b_Ndata_bb_hodo_tdc_te;   //!
   TBranch        *b_bb_hodo_tdc_te;   //!
   TBranch        *b_Ndata_bb_hodo_tdc_tot;   //!
   TBranch        *b_bb_hodo_tdc_tot;   //!
   TBranch        *b_Ndata_bb_hodo_tdccol;   //!
   TBranch        *b_bb_hodo_tdccol;   //!
   TBranch        *b_Ndata_bb_hodo_tdcelemID;   //!
   TBranch        *b_bb_hodo_tdcelemID;   //!
   TBranch        *b_Ndata_bb_hodo_tdclayer;   //!
   TBranch        *b_bb_hodo_tdclayer;   //!
   TBranch        *b_Ndata_bb_hodo_tdcrow;   //!
   TBranch        *b_bb_hodo_tdcrow;   //!
   TBranch        *b_Ndata_bb_prob_e;   //!
   TBranch        *b_bb_prob_e;   //!
   TBranch        *b_Ndata_bb_prob_pi;   //!
   TBranch        *b_bb_prob_pi;   //!
   TBranch        *b_Ndata_bb_ps_a;   //!
   TBranch        *b_bb_ps_a;   //!
   TBranch        *b_Ndata_bb_ps_a_amp;   //!
   TBranch        *b_bb_ps_a_amp;   //!
   TBranch        *b_Ndata_bb_ps_a_amp_p;   //!
   TBranch        *b_bb_ps_a_amp_p;   //!
   TBranch        *b_Ndata_bb_ps_a_c;   //!
   TBranch        *b_bb_ps_a_c;   //!
   TBranch        *b_Ndata_bb_ps_a_mult;   //!
   TBranch        *b_bb_ps_a_mult;   //!
   TBranch        *b_Ndata_bb_ps_a_p;   //!
   TBranch        *b_bb_ps_a_p;   //!
   TBranch        *b_Ndata_bb_ps_a_time;   //!
   TBranch        *b_bb_ps_a_time;   //!
   TBranch        *b_Ndata_bb_ps_adccol;   //!
   TBranch        *b_bb_ps_adccol;   //!
   TBranch        *b_Ndata_bb_ps_adcelemID;   //!
   TBranch        *b_bb_ps_adcelemID;   //!
   TBranch        *b_Ndata_bb_ps_adclayer;   //!
   TBranch        *b_bb_ps_adclayer;   //!
   TBranch        *b_Ndata_bb_ps_adcrow;   //!
   TBranch        *b_bb_ps_adcrow;   //!
   TBranch        *b_Ndata_bb_ps_clus_col;   //!
   TBranch        *b_bb_ps_clus_col;   //!
   TBranch        *b_Ndata_bb_ps_clus_e;   //!
   TBranch        *b_bb_ps_clus_e;   //!
   TBranch        *b_Ndata_bb_ps_clus_e_c;   //!
   TBranch        *b_bb_ps_clus_e_c;   //!
   TBranch        *b_Ndata_bb_ps_clus_eblk;   //!
   TBranch        *b_bb_ps_clus_eblk;   //!
   TBranch        *b_Ndata_bb_ps_clus_eblk_c;   //!
   TBranch        *b_bb_ps_clus_eblk_c;   //!
   TBranch        *b_Ndata_bb_ps_clus_id;   //!
   TBranch        *b_bb_ps_clus_id;   //!
   TBranch        *b_Ndata_bb_ps_clus_nblk;   //!
   TBranch        *b_bb_ps_clus_nblk;   //!
   TBranch        *b_Ndata_bb_ps_clus_row;   //!
   TBranch        *b_bb_ps_clus_row;   //!
   TBranch        *b_Ndata_bb_ps_clus_x;   //!
   TBranch        *b_bb_ps_clus_x;   //!
   TBranch        *b_Ndata_bb_ps_clus_y;   //!
   TBranch        *b_bb_ps_clus_y;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_col;   //!
   TBranch        *b_bb_ps_clus_blk_col;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_e;   //!
   TBranch        *b_bb_ps_clus_blk_e;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_e_c;   //!
   TBranch        *b_bb_ps_clus_blk_e_c;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_id;   //!
   TBranch        *b_bb_ps_clus_blk_id;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_row;   //!
   TBranch        *b_bb_ps_clus_blk_row;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_x;   //!
   TBranch        *b_bb_ps_clus_blk_x;   //!
   TBranch        *b_Ndata_bb_ps_clus_blk_y;   //!
   TBranch        *b_bb_ps_clus_blk_y;   //!
   TBranch        *b_Ndata_bb_ps_e_res;   //!
   TBranch        *b_bb_ps_e_res;   //!
   TBranch        *b_Ndata_bb_ps_goodblock_atime;   //!
   TBranch        *b_bb_ps_goodblock_atime;   //!
   TBranch        *b_Ndata_bb_ps_goodblock_col;   //!
   TBranch        *b_bb_ps_goodblock_col;   //!
   TBranch        *b_Ndata_bb_ps_goodblock_e;   //!
   TBranch        *b_bb_ps_goodblock_e;   //!
   TBranch        *b_Ndata_bb_ps_goodblock_id;   //!
   TBranch        *b_bb_ps_goodblock_id;   //!
   TBranch        *b_Ndata_bb_ps_goodblock_row;   //!
   TBranch        *b_bb_ps_goodblock_row;   //!
   TBranch        *b_Ndata_bb_ps_goodblock_x;   //!
   TBranch        *b_bb_ps_goodblock_x;   //!
   TBranch        *b_Ndata_bb_ps_goodblock_y;   //!
   TBranch        *b_bb_ps_goodblock_y;   //!
   TBranch        *b_Ndata_bb_ps_ped;   //!
   TBranch        *b_bb_ps_ped;   //!
   TBranch        *b_Ndata_bb_ps_x_res;   //!
   TBranch        *b_bb_ps_x_res;   //!
   TBranch        *b_Ndata_bb_ps_y_res;   //!
   TBranch        *b_bb_ps_y_res;   //!
   TBranch        *b_Ndata_bb_sh_a;   //!
   TBranch        *b_bb_sh_a;   //!
   TBranch        *b_Ndata_bb_sh_a_amp;   //!
   TBranch        *b_bb_sh_a_amp;   //!
   TBranch        *b_Ndata_bb_sh_a_amp_p;   //!
   TBranch        *b_bb_sh_a_amp_p;   //!
   TBranch        *b_Ndata_bb_sh_a_c;   //!
   TBranch        *b_bb_sh_a_c;   //!
   TBranch        *b_Ndata_bb_sh_a_mult;   //!
   TBranch        *b_bb_sh_a_mult;   //!
   TBranch        *b_Ndata_bb_sh_a_p;   //!
   TBranch        *b_bb_sh_a_p;   //!
   TBranch        *b_Ndata_bb_sh_a_time;   //!
   TBranch        *b_bb_sh_a_time;   //!
   TBranch        *b_Ndata_bb_sh_adccol;   //!
   TBranch        *b_bb_sh_adccol;   //!
   TBranch        *b_Ndata_bb_sh_adcelemID;   //!
   TBranch        *b_bb_sh_adcelemID;   //!
   TBranch        *b_Ndata_bb_sh_adclayer;   //!
   TBranch        *b_bb_sh_adclayer;   //!
   TBranch        *b_Ndata_bb_sh_adcrow;   //!
   TBranch        *b_bb_sh_adcrow;   //!
   TBranch        *b_Ndata_bb_sh_clus_col;   //!
   TBranch        *b_bb_sh_clus_col;   //!
   TBranch        *b_Ndata_bb_sh_clus_e;   //!
   TBranch        *b_bb_sh_clus_e;   //!
   TBranch        *b_Ndata_bb_sh_clus_e_c;   //!
   TBranch        *b_bb_sh_clus_e_c;   //!
   TBranch        *b_Ndata_bb_sh_clus_eblk;   //!
   TBranch        *b_bb_sh_clus_eblk;   //!
   TBranch        *b_Ndata_bb_sh_clus_eblk_c;   //!
   TBranch        *b_bb_sh_clus_eblk_c;   //!
   TBranch        *b_Ndata_bb_sh_clus_id;   //!
   TBranch        *b_bb_sh_clus_id;   //!
   TBranch        *b_Ndata_bb_sh_clus_nblk;   //!
   TBranch        *b_bb_sh_clus_nblk;   //!
   TBranch        *b_Ndata_bb_sh_clus_row;   //!
   TBranch        *b_bb_sh_clus_row;   //!
   TBranch        *b_Ndata_bb_sh_clus_x;   //!
   TBranch        *b_bb_sh_clus_x;   //!
   TBranch        *b_Ndata_bb_sh_clus_y;   //!
   TBranch        *b_bb_sh_clus_y;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_col;   //!
   TBranch        *b_bb_sh_clus_blk_col;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_e;   //!
   TBranch        *b_bb_sh_clus_blk_e;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_e_c;   //!
   TBranch        *b_bb_sh_clus_blk_e_c;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_id;   //!
   TBranch        *b_bb_sh_clus_blk_id;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_row;   //!
   TBranch        *b_bb_sh_clus_blk_row;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_x;   //!
   TBranch        *b_bb_sh_clus_blk_x;   //!
   TBranch        *b_Ndata_bb_sh_clus_blk_y;   //!
   TBranch        *b_bb_sh_clus_blk_y;   //!
   TBranch        *b_Ndata_bb_sh_e_res;   //!
   TBranch        *b_bb_sh_e_res;   //!
   TBranch        *b_Ndata_bb_sh_goodblock_atime;   //!
   TBranch        *b_bb_sh_goodblock_atime;   //!
   TBranch        *b_Ndata_bb_sh_goodblock_col;   //!
   TBranch        *b_bb_sh_goodblock_col;   //!
   TBranch        *b_Ndata_bb_sh_goodblock_e;   //!
   TBranch        *b_bb_sh_goodblock_e;   //!
   TBranch        *b_Ndata_bb_sh_goodblock_id;   //!
   TBranch        *b_bb_sh_goodblock_id;   //!
   TBranch        *b_Ndata_bb_sh_goodblock_row;   //!
   TBranch        *b_bb_sh_goodblock_row;   //!
   TBranch        *b_Ndata_bb_sh_goodblock_x;   //!
   TBranch        *b_bb_sh_goodblock_x;   //!
   TBranch        *b_Ndata_bb_sh_goodblock_y;   //!
   TBranch        *b_bb_sh_goodblock_y;   //!
   TBranch        *b_Ndata_bb_sh_ped;   //!
   TBranch        *b_bb_sh_ped;   //!
   TBranch        *b_Ndata_bb_sh_x_res;   //!
   TBranch        *b_bb_sh_x_res;   //!
   TBranch        *b_Ndata_bb_sh_y_res;   //!
   TBranch        *b_bb_sh_y_res;   //!
   TBranch        *b_Ndata_bb_tr_beta;   //!
   TBranch        *b_bb_tr_beta;   //!
   TBranch        *b_Ndata_bb_tr_chi2;   //!
   TBranch        *b_bb_tr_chi2;   //!
   TBranch        *b_Ndata_bb_tr_d_ph;   //!
   TBranch        *b_bb_tr_d_ph;   //!
   TBranch        *b_Ndata_bb_tr_d_th;   //!
   TBranch        *b_bb_tr_d_th;   //!
   TBranch        *b_Ndata_bb_tr_d_x;   //!
   TBranch        *b_bb_tr_d_x;   //!
   TBranch        *b_Ndata_bb_tr_d_y;   //!
   TBranch        *b_bb_tr_d_y;   //!
   TBranch        *b_Ndata_bb_tr_dbeta;   //!
   TBranch        *b_bb_tr_dbeta;   //!
   TBranch        *b_Ndata_bb_tr_dtime;   //!
   TBranch        *b_bb_tr_dtime;   //!
   TBranch        *b_Ndata_bb_tr_flag;   //!
   TBranch        *b_bb_tr_flag;   //!
   TBranch        *b_Ndata_bb_tr_ndof;   //!
   TBranch        *b_bb_tr_ndof;   //!
   TBranch        *b_Ndata_bb_tr_p;   //!
   TBranch        *b_bb_tr_p;   //!
   TBranch        *b_Ndata_bb_tr_pathl;   //!
   TBranch        *b_bb_tr_pathl;   //!
   TBranch        *b_Ndata_bb_tr_ph;   //!
   TBranch        *b_bb_tr_ph;   //!
   TBranch        *b_Ndata_bb_tr_px;   //!
   TBranch        *b_bb_tr_px;   //!
   TBranch        *b_Ndata_bb_tr_py;   //!
   TBranch        *b_bb_tr_py;   //!
   TBranch        *b_Ndata_bb_tr_pz;   //!
   TBranch        *b_bb_tr_pz;   //!
   TBranch        *b_Ndata_bb_tr_r_ph;   //!
   TBranch        *b_bb_tr_r_ph;   //!
   TBranch        *b_Ndata_bb_tr_r_th;   //!
   TBranch        *b_bb_tr_r_th;   //!
   TBranch        *b_Ndata_bb_tr_r_x;   //!
   TBranch        *b_bb_tr_r_x;   //!
   TBranch        *b_Ndata_bb_tr_r_y;   //!
   TBranch        *b_bb_tr_r_y;   //!
   TBranch        *b_Ndata_bb_tr_tg_dp;   //!
   TBranch        *b_bb_tr_tg_dp;   //!
   TBranch        *b_Ndata_bb_tr_tg_ph;   //!
   TBranch        *b_bb_tr_tg_ph;   //!
   TBranch        *b_Ndata_bb_tr_tg_th;   //!
   TBranch        *b_bb_tr_tg_th;   //!
   TBranch        *b_Ndata_bb_tr_tg_y;   //!
   TBranch        *b_bb_tr_tg_y;   //!
   TBranch        *b_Ndata_bb_tr_th;   //!
   TBranch        *b_bb_tr_th;   //!
   TBranch        *b_Ndata_bb_tr_time;   //!
   TBranch        *b_bb_tr_time;   //!
   TBranch        *b_Ndata_bb_tr_vx;   //!
   TBranch        *b_bb_tr_vx;   //!
   TBranch        *b_Ndata_bb_tr_vy;   //!
   TBranch        *b_bb_tr_vy;   //!
   TBranch        *b_Ndata_bb_tr_vz;   //!
   TBranch        *b_bb_tr_vz;   //!
   TBranch        *b_Ndata_bb_tr_x;   //!
   TBranch        *b_bb_tr_x;   //!
   TBranch        *b_Ndata_bb_tr_y;   //!
   TBranch        *b_bb_tr_y;   //!
   TBranch        *b_Ndata_bb_x_bcp;   //!
   TBranch        *b_bb_x_bcp;   //!
   TBranch        *b_Ndata_bb_x_fcp;   //!
   TBranch        *b_bb_x_fcp;   //!
   TBranch        *b_Ndata_bb_y_bcp;   //!
   TBranch        *b_bb_y_bcp;   //!
   TBranch        *b_Ndata_bb_y_fcp;   //!
   TBranch        *b_bb_y_fcp;   //!
   TBranch        *b_Ndata_sbs_hcal_a;   //!
   TBranch        *b_sbs_hcal_a;   //!
   TBranch        *b_Ndata_sbs_hcal_a_amp;   //!
   TBranch        *b_sbs_hcal_a_amp;   //!
   TBranch        *b_Ndata_sbs_hcal_a_amp_p;   //!
   TBranch        *b_sbs_hcal_a_amp_p;   //!
   TBranch        *b_Ndata_sbs_hcal_a_c;   //!
   TBranch        *b_sbs_hcal_a_c;   //!
   TBranch        *b_Ndata_sbs_hcal_a_mult;   //!
   TBranch        *b_sbs_hcal_a_mult;   //!
   TBranch        *b_Ndata_sbs_hcal_a_p;   //!
   TBranch        *b_sbs_hcal_a_p;   //!
   TBranch        *b_Ndata_sbs_hcal_a_time;   //!
   TBranch        *b_sbs_hcal_a_time;   //!
   TBranch        *b_Ndata_sbs_hcal_adccol;   //!
   TBranch        *b_sbs_hcal_adccol;   //!
   TBranch        *b_Ndata_sbs_hcal_adcelemID;   //!
   TBranch        *b_sbs_hcal_adcelemID;   //!
   TBranch        *b_Ndata_sbs_hcal_adclayer;   //!
   TBranch        *b_sbs_hcal_adclayer;   //!
   TBranch        *b_Ndata_sbs_hcal_adcrow;   //!
   TBranch        *b_sbs_hcal_adcrow;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_col;   //!
   TBranch        *b_sbs_hcal_clus_col;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_e;   //!
   TBranch        *b_sbs_hcal_clus_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_e_c;   //!
   TBranch        *b_sbs_hcal_clus_e_c;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_eblk;   //!
   TBranch        *b_sbs_hcal_clus_eblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_eblk_c;   //!
   TBranch        *b_sbs_hcal_clus_eblk_c;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_id;   //!
   TBranch        *b_sbs_hcal_clus_id;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_nblk;   //!
   TBranch        *b_sbs_hcal_clus_nblk;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_row;   //!
   TBranch        *b_sbs_hcal_clus_row;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_x;   //!
   TBranch        *b_sbs_hcal_clus_x;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_y;   //!
   TBranch        *b_sbs_hcal_clus_y;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_col;   //!
   TBranch        *b_sbs_hcal_clus_blk_col;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_e;   //!
   TBranch        *b_sbs_hcal_clus_blk_e;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_e_c;   //!
   TBranch        *b_sbs_hcal_clus_blk_e_c;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_id;   //!
   TBranch        *b_sbs_hcal_clus_blk_id;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_row;   //!
   TBranch        *b_sbs_hcal_clus_blk_row;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_x;   //!
   TBranch        *b_sbs_hcal_clus_blk_x;   //!
   TBranch        *b_Ndata_sbs_hcal_clus_blk_y;   //!
   TBranch        *b_sbs_hcal_clus_blk_y;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_atime;   //!
   TBranch        *b_sbs_hcal_goodblock_atime;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_col;   //!
   TBranch        *b_sbs_hcal_goodblock_col;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_e;   //!
   TBranch        *b_sbs_hcal_goodblock_e;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_id;   //!
   TBranch        *b_sbs_hcal_goodblock_id;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_row;   //!
   TBranch        *b_sbs_hcal_goodblock_row;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_x;   //!
   TBranch        *b_sbs_hcal_goodblock_x;   //!
   TBranch        *b_Ndata_sbs_hcal_goodblock_y;   //!
   TBranch        *b_sbs_hcal_goodblock_y;   //!
   TBranch        *b_Ndata_sbs_hcal_ped;   //!
   TBranch        *b_sbs_hcal_ped;   //!
   TBranch        *b_Ndata_sbs_hcal_tdc;   //!
   TBranch        *b_sbs_hcal_tdc;   //!
   TBranch        *b_Ndata_sbs_hcal_tdc_mult;   //!
   TBranch        *b_sbs_hcal_tdc_mult;   //!
   TBranch        *b_Ndata_sbs_hcal_tdccol;   //!
   TBranch        *b_sbs_hcal_tdccol;   //!
   TBranch        *b_Ndata_sbs_hcal_tdcelemID;   //!
   TBranch        *b_sbs_hcal_tdcelemID;   //!
   TBranch        *b_Ndata_sbs_hcal_tdclayer;   //!
   TBranch        *b_sbs_hcal_tdclayer;   //!
   TBranch        *b_Ndata_sbs_hcal_tdcrow;   //!
   TBranch        *b_sbs_hcal_tdcrow;   //!
   TBranch        *b_Ndata_sbs_tr_beta;   //!
   TBranch        *b_sbs_tr_beta;   //!
   TBranch        *b_Ndata_sbs_tr_chi2;   //!
   TBranch        *b_sbs_tr_chi2;   //!
   TBranch        *b_Ndata_sbs_tr_d_ph;   //!
   TBranch        *b_sbs_tr_d_ph;   //!
   TBranch        *b_Ndata_sbs_tr_d_th;   //!
   TBranch        *b_sbs_tr_d_th;   //!
   TBranch        *b_Ndata_sbs_tr_d_x;   //!
   TBranch        *b_sbs_tr_d_x;   //!
   TBranch        *b_Ndata_sbs_tr_d_y;   //!
   TBranch        *b_sbs_tr_d_y;   //!
   TBranch        *b_Ndata_sbs_tr_dbeta;   //!
   TBranch        *b_sbs_tr_dbeta;   //!
   TBranch        *b_Ndata_sbs_tr_dtime;   //!
   TBranch        *b_sbs_tr_dtime;   //!
   TBranch        *b_Ndata_sbs_tr_flag;   //!
   TBranch        *b_sbs_tr_flag;   //!
   TBranch        *b_Ndata_sbs_tr_ndof;   //!
   TBranch        *b_sbs_tr_ndof;   //!
   TBranch        *b_Ndata_sbs_tr_p;   //!
   TBranch        *b_sbs_tr_p;   //!
   TBranch        *b_Ndata_sbs_tr_pathl;   //!
   TBranch        *b_sbs_tr_pathl;   //!
   TBranch        *b_Ndata_sbs_tr_ph;   //!
   TBranch        *b_sbs_tr_ph;   //!
   TBranch        *b_Ndata_sbs_tr_px;   //!
   TBranch        *b_sbs_tr_px;   //!
   TBranch        *b_Ndata_sbs_tr_py;   //!
   TBranch        *b_sbs_tr_py;   //!
   TBranch        *b_Ndata_sbs_tr_pz;   //!
   TBranch        *b_sbs_tr_pz;   //!
   TBranch        *b_Ndata_sbs_tr_r_ph;   //!
   TBranch        *b_sbs_tr_r_ph;   //!
   TBranch        *b_Ndata_sbs_tr_r_th;   //!
   TBranch        *b_sbs_tr_r_th;   //!
   TBranch        *b_Ndata_sbs_tr_r_x;   //!
   TBranch        *b_sbs_tr_r_x;   //!
   TBranch        *b_Ndata_sbs_tr_r_y;   //!
   TBranch        *b_sbs_tr_r_y;   //!
   TBranch        *b_Ndata_sbs_tr_tg_dp;   //!
   TBranch        *b_sbs_tr_tg_dp;   //!
   TBranch        *b_Ndata_sbs_tr_tg_ph;   //!
   TBranch        *b_sbs_tr_tg_ph;   //!
   TBranch        *b_Ndata_sbs_tr_tg_th;   //!
   TBranch        *b_sbs_tr_tg_th;   //!
   TBranch        *b_Ndata_sbs_tr_tg_y;   //!
   TBranch        *b_sbs_tr_tg_y;   //!
   TBranch        *b_Ndata_sbs_tr_th;   //!
   TBranch        *b_sbs_tr_th;   //!
   TBranch        *b_Ndata_sbs_tr_time;   //!
   TBranch        *b_sbs_tr_time;   //!
   TBranch        *b_Ndata_sbs_tr_vx;   //!
   TBranch        *b_sbs_tr_vx;   //!
   TBranch        *b_Ndata_sbs_tr_vy;   //!
   TBranch        *b_sbs_tr_vy;   //!
   TBranch        *b_Ndata_sbs_tr_vz;   //!
   TBranch        *b_sbs_tr_vz;   //!
   TBranch        *b_Ndata_sbs_tr_x;   //!
   TBranch        *b_sbs_tr_x;   //!
   TBranch        *b_Ndata_sbs_tr_y;   //!
   TBranch        *b_sbs_tr_y;   //!
   TBranch        *b_MC_hit_n;   //!
   TBranch        *b_MC_mc_px;   //!
   TBranch        *b_MC_mc_py;   //!
   TBranch        *b_MC_mc_pz;   //!
   TBranch        *b_MC_mc_vx;   //!
   TBranch        *b_MC_mc_vy;   //!
   TBranch        *b_MC_mc_vz;   //!
   TBranch        *b_MC_nbbgemhits;   //!
   TBranch        *b_MC_nbbtracks;   //!
   TBranch        *b_MC_pt_n;   //!
   TBranch        *b_MC_tr_n;   //!
   TBranch        *b_MC_weight;   //!
   TBranch        *b_bb_gem_hit_ngoodhits;   //!
   TBranch        *b_bb_gem_nlayershit;   //!
   TBranch        *b_bb_gem_nlayershitu;   //!
   TBranch        *b_bb_gem_nlayershituv;   //!
   TBranch        *b_bb_gem_nlayershitv;   //!
   TBranch        *b_bb_gem_track_besttrack;   //!
   TBranch        *b_bb_gem_track_ntrack;   //!
   TBranch        *b_bb_grinch_nclus;   //!
   TBranch        *b_bb_grinch_nhits;   //!
   TBranch        *b_bb_hodo_bar_ngoodbars;   //!
   TBranch        *b_bb_hodo_nclus;   //!
   TBranch        *b_bb_hodo_ngoodADChits;   //!
   TBranch        *b_bb_hodo_ngoodTDChits;   //!
   TBranch        *b_bb_hodo_nhits;   //!
   TBranch        *b_bb_hodo_nrefhits;   //!
   TBranch        *b_bb_ps_colblk;   //!
   TBranch        *b_bb_ps_e;   //!
   TBranch        *b_bb_ps_e_c;   //!
   TBranch        *b_bb_ps_e_m_res;   //!
   TBranch        *b_bb_ps_eblk;   //!
   TBranch        *b_bb_ps_eblk_c;   //!
   TBranch        *b_bb_ps_idblk;   //!
   TBranch        *b_bb_ps_nblk;   //!
   TBranch        *b_bb_ps_nclus;   //!
   TBranch        *b_bb_ps_ngoodADChits;   //!
   TBranch        *b_bb_ps_ngoodTDChits;   //!
   TBranch        *b_bb_ps_nhits;   //!
   TBranch        *b_bb_ps_nrefhits;   //!
   TBranch        *b_bb_ps_rowblk;   //!
   TBranch        *b_bb_ps_x;   //!
   TBranch        *b_bb_ps_x_m_res;   //!
   TBranch        *b_bb_ps_y;   //!
   TBranch        *b_bb_ps_y_m_res;   //!
   TBranch        *b_bb_sh_colblk;   //!
   TBranch        *b_bb_sh_e;   //!
   TBranch        *b_bb_sh_e_c;   //!
   TBranch        *b_bb_sh_e_m_res;   //!
   TBranch        *b_bb_sh_eblk;   //!
   TBranch        *b_bb_sh_eblk_c;   //!
   TBranch        *b_bb_sh_idblk;   //!
   TBranch        *b_bb_sh_nblk;   //!
   TBranch        *b_bb_sh_nclus;   //!
   TBranch        *b_bb_sh_ngoodADChits;   //!
   TBranch        *b_bb_sh_ngoodTDChits;   //!
   TBranch        *b_bb_sh_nhits;   //!
   TBranch        *b_bb_sh_nrefhits;   //!
   TBranch        *b_bb_sh_rowblk;   //!
   TBranch        *b_bb_sh_x;   //!
   TBranch        *b_bb_sh_x_m_res;   //!
   TBranch        *b_bb_sh_y;   //!
   TBranch        *b_bb_sh_y_m_res;   //!
   TBranch        *b_bb_tr_n;   //!
   TBranch        *b_sbs_hcal_colblk;   //!
   TBranch        *b_sbs_hcal_e;   //!
   TBranch        *b_sbs_hcal_e_c;   //!
   TBranch        *b_sbs_hcal_eblk;   //!
   TBranch        *b_sbs_hcal_eblk_c;   //!
   TBranch        *b_sbs_hcal_idblk;   //!
   TBranch        *b_sbs_hcal_ledbit;   //!
   TBranch        *b_sbs_hcal_ledcount;   //!
   TBranch        *b_sbs_hcal_nblk;   //!
   TBranch        *b_sbs_hcal_nclus;   //!
   TBranch        *b_sbs_hcal_ngoodADChits;   //!
   TBranch        *b_sbs_hcal_ngoodTDChits;   //!
   TBranch        *b_sbs_hcal_nhits;   //!
   TBranch        *b_sbs_hcal_nrefhits;   //!
   TBranch        *b_sbs_hcal_rowblk;   //!
   TBranch        *b_sbs_hcal_x;   //!
   TBranch        *b_sbs_hcal_y;   //!
   TBranch        *b_sbs_status;   //!
   TBranch        *b_sbs_tr_n;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtTime;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtNum;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtType;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fEvtLen;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fHelicity;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fTargetPol;   //!
   TBranch        *b_Event_Branch_fEvtHdr_fRun;   //!

   gmn_rec_tree(TTree *tree=0);
   virtual ~gmn_rec_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gmn_rec_tree_cxx
gmn_rec_tree::gmn_rec_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("replayed_simdigtest_2_20211004.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("replayed_simdigtest_2_20211004.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

gmn_rec_tree::~gmn_rec_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gmn_rec_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gmn_rec_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gmn_rec_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Ndata.MC.bbgemhit_edep", &Ndata_MC_bbgemhit_edep, &b_Ndata_MC_bbgemhit_edep);
   fChain->SetBranchAddress("MC.bbgemhit_edep", MC_bbgemhit_edep, &b_MC_bbgemhit_edep);
   fChain->SetBranchAddress("Ndata.MC.bbgemhit_mid", &Ndata_MC_bbgemhit_mid, &b_Ndata_MC_bbgemhit_mid);
   fChain->SetBranchAddress("MC.bbgemhit_mid", MC_bbgemhit_mid, &b_MC_bbgemhit_mid);
   fChain->SetBranchAddress("Ndata.MC.bbgemhit_pid", &Ndata_MC_bbgemhit_pid, &b_Ndata_MC_bbgemhit_pid);
   fChain->SetBranchAddress("MC.bbgemhit_pid", MC_bbgemhit_pid, &b_MC_bbgemhit_pid);
   fChain->SetBranchAddress("Ndata.MC.bbgemhit_plane", &Ndata_MC_bbgemhit_plane, &b_Ndata_MC_bbgemhit_plane);
   fChain->SetBranchAddress("MC.bbgemhit_plane", MC_bbgemhit_plane, &b_MC_bbgemhit_plane);
   fChain->SetBranchAddress("Ndata.MC.bbgemhit_tid", &Ndata_MC_bbgemhit_tid, &b_Ndata_MC_bbgemhit_tid);
   fChain->SetBranchAddress("MC.bbgemhit_tid", MC_bbgemhit_tid, &b_MC_bbgemhit_tid);
   fChain->SetBranchAddress("Ndata.MC.bbgemhit_x", &Ndata_MC_bbgemhit_x, &b_Ndata_MC_bbgemhit_x);
   fChain->SetBranchAddress("MC.bbgemhit_x", MC_bbgemhit_x, &b_MC_bbgemhit_x);
   fChain->SetBranchAddress("Ndata.MC.bbgemhit_y", &Ndata_MC_bbgemhit_y, &b_Ndata_MC_bbgemhit_y);
   fChain->SetBranchAddress("MC.bbgemhit_y", MC_bbgemhit_y, &b_MC_bbgemhit_y);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_dx", &Ndata_MC_bbtrack_dx, &b_Ndata_MC_bbtrack_dx);
   fChain->SetBranchAddress("MC.bbtrack_dx", MC_bbtrack_dx, &b_MC_bbtrack_dx);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_dy", &Ndata_MC_bbtrack_dy, &b_Ndata_MC_bbtrack_dy);
   fChain->SetBranchAddress("MC.bbtrack_dy", MC_bbtrack_dy, &b_MC_bbtrack_dy);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_mid", &Ndata_MC_bbtrack_mid, &b_Ndata_MC_bbtrack_mid);
   fChain->SetBranchAddress("MC.bbtrack_mid", MC_bbtrack_mid, &b_MC_bbtrack_mid);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_nhits", &Ndata_MC_bbtrack_nhits, &b_Ndata_MC_bbtrack_nhits);
   fChain->SetBranchAddress("MC.bbtrack_nhits", MC_bbtrack_nhits, &b_MC_bbtrack_nhits);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_p", &Ndata_MC_bbtrack_p, &b_Ndata_MC_bbtrack_p);
   fChain->SetBranchAddress("MC.bbtrack_p", MC_bbtrack_p, &b_MC_bbtrack_p);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_pid", &Ndata_MC_bbtrack_pid, &b_Ndata_MC_bbtrack_pid);
   fChain->SetBranchAddress("MC.bbtrack_pid", MC_bbtrack_pid, &b_MC_bbtrack_pid);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_tid", &Ndata_MC_bbtrack_tid, &b_Ndata_MC_bbtrack_tid);
   fChain->SetBranchAddress("MC.bbtrack_tid", MC_bbtrack_tid, &b_MC_bbtrack_tid);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_x", &Ndata_MC_bbtrack_x, &b_Ndata_MC_bbtrack_x);
   fChain->SetBranchAddress("MC.bbtrack_x", MC_bbtrack_x, &b_MC_bbtrack_x);
   fChain->SetBranchAddress("Ndata.MC.bbtrack_y", &Ndata_MC_bbtrack_y, &b_Ndata_MC_bbtrack_y);
   fChain->SetBranchAddress("MC.bbtrack_y", MC_bbtrack_y, &b_MC_bbtrack_y);
   fChain->SetBranchAddress("Ndata.MC.pt.clustsz", &Ndata_MC_pt_clustsz, &b_Ndata_MC_pt_clustsz);
   fChain->SetBranchAddress("MC.pt.clustsz", &MC_pt_clustsz, &b_MC_pt_clustsz);
   fChain->SetBranchAddress("Ndata.MC.pt.deflect", &Ndata_MC_pt_deflect, &b_Ndata_MC_pt_deflect);
   fChain->SetBranchAddress("MC.pt.deflect", &MC_pt_deflect, &b_MC_pt_deflect);
   fChain->SetBranchAddress("Ndata.MC.pt.deltaE", &Ndata_MC_pt_deltaE, &b_Ndata_MC_pt_deltaE);
   fChain->SetBranchAddress("MC.pt.deltaE", &MC_pt_deltaE, &b_MC_pt_deltaE);
   fChain->SetBranchAddress("Ndata.MC.pt.hitres", &Ndata_MC_pt_hitres, &b_Ndata_MC_pt_hitres);
   fChain->SetBranchAddress("MC.pt.hitres", &MC_pt_hitres, &b_MC_pt_hitres);
   fChain->SetBranchAddress("Ndata.MC.pt.nfound", &Ndata_MC_pt_nfound, &b_Ndata_MC_pt_nfound);
   fChain->SetBranchAddress("MC.pt.nfound", &MC_pt_nfound, &b_MC_pt_nfound);
   fChain->SetBranchAddress("Ndata.MC.pt.p", &Ndata_MC_pt_p, &b_Ndata_MC_pt_p);
   fChain->SetBranchAddress("MC.pt.p", &MC_pt_p, &b_MC_pt_p);
   fChain->SetBranchAddress("Ndata.MC.pt.ph", &Ndata_MC_pt_ph, &b_Ndata_MC_pt_ph);
   fChain->SetBranchAddress("MC.pt.ph", &MC_pt_ph, &b_MC_pt_ph);
   fChain->SetBranchAddress("Ndata.MC.pt.phdir", &Ndata_MC_pt_phdir, &b_Ndata_MC_pt_phdir);
   fChain->SetBranchAddress("MC.pt.phdir", &MC_pt_phdir, &b_MC_pt_phdir);
   fChain->SetBranchAddress("Ndata.MC.pt.phi", &Ndata_MC_pt_phi, &b_Ndata_MC_pt_phi);
   fChain->SetBranchAddress("MC.pt.phi", &MC_pt_phi, &b_MC_pt_phi);
   fChain->SetBranchAddress("Ndata.MC.pt.plane", &Ndata_MC_pt_plane, &b_Ndata_MC_pt_plane);
   fChain->SetBranchAddress("MC.pt.plane", &MC_pt_plane, &b_MC_pt_plane);
   fChain->SetBranchAddress("Ndata.MC.pt.r", &Ndata_MC_pt_r, &b_Ndata_MC_pt_r);
   fChain->SetBranchAddress("MC.pt.r", &MC_pt_r, &b_MC_pt_r);
   fChain->SetBranchAddress("Ndata.MC.pt.status", &Ndata_MC_pt_status, &b_Ndata_MC_pt_status);
   fChain->SetBranchAddress("MC.pt.status", &MC_pt_status, &b_MC_pt_status);
   fChain->SetBranchAddress("Ndata.MC.pt.th", &Ndata_MC_pt_th, &b_Ndata_MC_pt_th);
   fChain->SetBranchAddress("MC.pt.th", &MC_pt_th, &b_MC_pt_th);
   fChain->SetBranchAddress("Ndata.MC.pt.thdir", &Ndata_MC_pt_thdir, &b_Ndata_MC_pt_thdir);
   fChain->SetBranchAddress("MC.pt.thdir", &MC_pt_thdir, &b_MC_pt_thdir);
   fChain->SetBranchAddress("Ndata.MC.pt.theta", &Ndata_MC_pt_theta, &b_Ndata_MC_pt_theta);
   fChain->SetBranchAddress("MC.pt.theta", &MC_pt_theta, &b_MC_pt_theta);
   fChain->SetBranchAddress("Ndata.MC.pt.time", &Ndata_MC_pt_time, &b_Ndata_MC_pt_time);
   fChain->SetBranchAddress("MC.pt.time", &MC_pt_time, &b_MC_pt_time);
   fChain->SetBranchAddress("Ndata.MC.pt.tof", &Ndata_MC_pt_tof, &b_Ndata_MC_pt_tof);
   fChain->SetBranchAddress("MC.pt.tof", &MC_pt_tof, &b_MC_pt_tof);
   fChain->SetBranchAddress("Ndata.MC.pt.trkres", &Ndata_MC_pt_trkres, &b_Ndata_MC_pt_trkres);
   fChain->SetBranchAddress("MC.pt.trkres", &MC_pt_trkres, &b_MC_pt_trkres);
   fChain->SetBranchAddress("Ndata.MC.pt.type", &Ndata_MC_pt_type, &b_Ndata_MC_pt_type);
   fChain->SetBranchAddress("MC.pt.type", &MC_pt_type, &b_MC_pt_type);
   fChain->SetBranchAddress("Ndata.MC.pt.x", &Ndata_MC_pt_x, &b_Ndata_MC_pt_x);
   fChain->SetBranchAddress("MC.pt.x", &MC_pt_x, &b_MC_pt_x);
   fChain->SetBranchAddress("Ndata.MC.pt.y", &Ndata_MC_pt_y, &b_Ndata_MC_pt_y);
   fChain->SetBranchAddress("MC.pt.y", &MC_pt_y, &b_MC_pt_y);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCU", &Ndata_bb_gem_hit_ADCU, &b_Ndata_bb_gem_hit_ADCU);
   fChain->SetBranchAddress("bb.gem.hit.ADCU", bb_gem_hit_ADCU, &b_bb_gem_hit_ADCU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCV", &Ndata_bb_gem_hit_ADCV, &b_Ndata_bb_gem_hit_ADCV);
   fChain->SetBranchAddress("bb.gem.hit.ADCV", bb_gem_hit_ADCV, &b_bb_gem_hit_ADCV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCasym", &Ndata_bb_gem_hit_ADCasym, &b_Ndata_bb_gem_hit_ADCasym);
   fChain->SetBranchAddress("bb.gem.hit.ADCasym", bb_gem_hit_ADCasym, &b_bb_gem_hit_ADCasym);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCavg", &Ndata_bb_gem_hit_ADCavg, &b_Ndata_bb_gem_hit_ADCavg);
   fChain->SetBranchAddress("bb.gem.hit.ADCavg", bb_gem_hit_ADCavg, &b_bb_gem_hit_ADCavg);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampU", &Ndata_bb_gem_hit_ADCmaxsampU, &b_Ndata_bb_gem_hit_ADCmaxsampU);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampU", bb_gem_hit_ADCmaxsampU, &b_bb_gem_hit_ADCmaxsampU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampUclust", &Ndata_bb_gem_hit_ADCmaxsampUclust, &b_Ndata_bb_gem_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampUclust", bb_gem_hit_ADCmaxsampUclust, &b_bb_gem_hit_ADCmaxsampUclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampV", &Ndata_bb_gem_hit_ADCmaxsampV, &b_Ndata_bb_gem_hit_ADCmaxsampV);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampV", bb_gem_hit_ADCmaxsampV, &b_bb_gem_hit_ADCmaxsampV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxsampVclust", &Ndata_bb_gem_hit_ADCmaxsampVclust, &b_Ndata_bb_gem_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxsampVclust", bb_gem_hit_ADCmaxsampVclust, &b_bb_gem_hit_ADCmaxsampVclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxstripU", &Ndata_bb_gem_hit_ADCmaxstripU, &b_Ndata_bb_gem_hit_ADCmaxstripU);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxstripU", bb_gem_hit_ADCmaxstripU, &b_bb_gem_hit_ADCmaxstripU);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ADCmaxstripV", &Ndata_bb_gem_hit_ADCmaxstripV, &b_Ndata_bb_gem_hit_ADCmaxstripV);
   fChain->SetBranchAddress("bb.gem.hit.ADCmaxstripV", bb_gem_hit_ADCmaxstripV, &b_bb_gem_hit_ADCmaxstripV);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_U", &Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U, &b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("bb.gem.hit.BUILD_ALL_SAMPLES_U", bb_gem_hit_BUILD_ALL_SAMPLES_U, &b_bb_gem_hit_BUILD_ALL_SAMPLES_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_V", &Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V, &b_Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("bb.gem.hit.BUILD_ALL_SAMPLES_V", bb_gem_hit_BUILD_ALL_SAMPLES_V, &b_bb_gem_hit_BUILD_ALL_SAMPLES_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.CM_OR_U", &Ndata_bb_gem_hit_CM_OR_U, &b_Ndata_bb_gem_hit_CM_OR_U);
   fChain->SetBranchAddress("bb.gem.hit.CM_OR_U", bb_gem_hit_CM_OR_U, &b_bb_gem_hit_CM_OR_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.CM_OR_V", &Ndata_bb_gem_hit_CM_OR_V, &b_Ndata_bb_gem_hit_CM_OR_V);
   fChain->SetBranchAddress("bb.gem.hit.CM_OR_V", bb_gem_hit_CM_OR_V, &b_bb_gem_hit_CM_OR_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ENABLE_CM_U", &Ndata_bb_gem_hit_ENABLE_CM_U, &b_Ndata_bb_gem_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("bb.gem.hit.ENABLE_CM_U", bb_gem_hit_ENABLE_CM_U, &b_bb_gem_hit_ENABLE_CM_U);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ENABLE_CM_V", &Ndata_bb_gem_hit_ENABLE_CM_V, &b_Ndata_bb_gem_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("bb.gem.hit.ENABLE_CM_V", bb_gem_hit_ENABLE_CM_V, &b_bb_gem_hit_ENABLE_CM_V);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Tavg", &Ndata_bb_gem_hit_Tavg, &b_Ndata_bb_gem_hit_Tavg);
   fChain->SetBranchAddress("bb.gem.hit.Tavg", bb_gem_hit_Tavg, &b_bb_gem_hit_Tavg);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Utime", &Ndata_bb_gem_hit_Utime, &b_Ndata_bb_gem_hit_Utime);
   fChain->SetBranchAddress("bb.gem.hit.Utime", bb_gem_hit_Utime, &b_bb_gem_hit_Utime);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.UtimeMaxStrip", &Ndata_bb_gem_hit_UtimeMaxStrip, &b_Ndata_bb_gem_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("bb.gem.hit.UtimeMaxStrip", bb_gem_hit_UtimeMaxStrip, &b_bb_gem_hit_UtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.Vtime", &Ndata_bb_gem_hit_Vtime, &b_Ndata_bb_gem_hit_Vtime);
   fChain->SetBranchAddress("bb.gem.hit.Vtime", bb_gem_hit_Vtime, &b_bb_gem_hit_Vtime);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.VtimeMaxStrip", &Ndata_bb_gem_hit_VtimeMaxStrip, &b_Ndata_bb_gem_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("bb.gem.hit.VtimeMaxStrip", bb_gem_hit_VtimeMaxStrip, &b_bb_gem_hit_VtimeMaxStrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ccor_clust", &Ndata_bb_gem_hit_ccor_clust, &b_Ndata_bb_gem_hit_ccor_clust);
   fChain->SetBranchAddress("bb.gem.hit.ccor_clust", bb_gem_hit_ccor_clust, &b_bb_gem_hit_ccor_clust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ccor_strip", &Ndata_bb_gem_hit_ccor_strip, &b_Ndata_bb_gem_hit_ccor_strip);
   fChain->SetBranchAddress("bb.gem.hit.ccor_strip", bb_gem_hit_ccor_strip, &b_bb_gem_hit_ccor_strip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.deltat", &Ndata_bb_gem_hit_deltat, &b_Ndata_bb_gem_hit_deltat);
   fChain->SetBranchAddress("bb.gem.hit.deltat", bb_gem_hit_deltat, &b_bb_gem_hit_deltat);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.eresidu", &Ndata_bb_gem_hit_eresidu, &b_Ndata_bb_gem_hit_eresidu);
   fChain->SetBranchAddress("bb.gem.hit.eresidu", bb_gem_hit_eresidu, &b_bb_gem_hit_eresidu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.eresidv", &Ndata_bb_gem_hit_eresidv, &b_Ndata_bb_gem_hit_eresidv);
   fChain->SetBranchAddress("bb.gem.hit.eresidv", bb_gem_hit_eresidv, &b_bb_gem_hit_eresidv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUclust", &Ndata_bb_gem_hit_isampmaxUclust, &b_Ndata_bb_gem_hit_isampmaxUclust);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxUclust", bb_gem_hit_isampmaxUclust, &b_bb_gem_hit_isampmaxUclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxUstrip", &Ndata_bb_gem_hit_isampmaxUstrip, &b_Ndata_bb_gem_hit_isampmaxUstrip);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxUstrip", bb_gem_hit_isampmaxUstrip, &b_bb_gem_hit_isampmaxUstrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxVclust", &Ndata_bb_gem_hit_isampmaxVclust, &b_Ndata_bb_gem_hit_isampmaxVclust);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxVclust", bb_gem_hit_isampmaxVclust, &b_bb_gem_hit_isampmaxVclust);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.isampmaxVstrip", &Ndata_bb_gem_hit_isampmaxVstrip, &b_Ndata_bb_gem_hit_isampmaxVstrip);
   fChain->SetBranchAddress("bb.gem.hit.isampmaxVstrip", bb_gem_hit_isampmaxVstrip, &b_bb_gem_hit_isampmaxVstrip);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.layer", &Ndata_bb_gem_hit_layer, &b_Ndata_bb_gem_hit_layer);
   fChain->SetBranchAddress("bb.gem.hit.layer", bb_gem_hit_layer, &b_bb_gem_hit_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.module", &Ndata_bb_gem_hit_module, &b_Ndata_bb_gem_hit_module);
   fChain->SetBranchAddress("bb.gem.hit.module", bb_gem_hit_module, &b_bb_gem_hit_module);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.nstripu", &Ndata_bb_gem_hit_nstripu, &b_Ndata_bb_gem_hit_nstripu);
   fChain->SetBranchAddress("bb.gem.hit.nstripu", bb_gem_hit_nstripu, &b_bb_gem_hit_nstripu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.nstripv", &Ndata_bb_gem_hit_nstripv, &b_Ndata_bb_gem_hit_nstripv);
   fChain->SetBranchAddress("bb.gem.hit.nstripv", bb_gem_hit_nstripv, &b_bb_gem_hit_nstripv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.residu", &Ndata_bb_gem_hit_residu, &b_Ndata_bb_gem_hit_residu);
   fChain->SetBranchAddress("bb.gem.hit.residu", bb_gem_hit_residu, &b_bb_gem_hit_residu);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.residv", &Ndata_bb_gem_hit_residv, &b_Ndata_bb_gem_hit_residv);
   fChain->SetBranchAddress("bb.gem.hit.residv", bb_gem_hit_residv, &b_bb_gem_hit_residv);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.trackindex", &Ndata_bb_gem_hit_trackindex, &b_Ndata_bb_gem_hit_trackindex);
   fChain->SetBranchAddress("bb.gem.hit.trackindex", bb_gem_hit_trackindex, &b_bb_gem_hit_trackindex);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.u", &Ndata_bb_gem_hit_u, &b_Ndata_bb_gem_hit_u);
   fChain->SetBranchAddress("bb.gem.hit.u", bb_gem_hit_u, &b_bb_gem_hit_u);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.umoment", &Ndata_bb_gem_hit_umoment, &b_Ndata_bb_gem_hit_umoment);
   fChain->SetBranchAddress("bb.gem.hit.umoment", bb_gem_hit_umoment, &b_bb_gem_hit_umoment);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.usigma", &Ndata_bb_gem_hit_usigma, &b_Ndata_bb_gem_hit_usigma);
   fChain->SetBranchAddress("bb.gem.hit.usigma", bb_gem_hit_usigma, &b_bb_gem_hit_usigma);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustriphi", &Ndata_bb_gem_hit_ustriphi, &b_Ndata_bb_gem_hit_ustriphi);
   fChain->SetBranchAddress("bb.gem.hit.ustriphi", bb_gem_hit_ustriphi, &b_bb_gem_hit_ustriphi);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustriplo", &Ndata_bb_gem_hit_ustriplo, &b_Ndata_bb_gem_hit_ustriplo);
   fChain->SetBranchAddress("bb.gem.hit.ustriplo", bb_gem_hit_ustriplo, &b_bb_gem_hit_ustriplo);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ustripmax", &Ndata_bb_gem_hit_ustripmax, &b_Ndata_bb_gem_hit_ustripmax);
   fChain->SetBranchAddress("bb.gem.hit.ustripmax", bb_gem_hit_ustripmax, &b_bb_gem_hit_ustripmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.v", &Ndata_bb_gem_hit_v, &b_Ndata_bb_gem_hit_v);
   fChain->SetBranchAddress("bb.gem.hit.v", bb_gem_hit_v, &b_bb_gem_hit_v);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vmoment", &Ndata_bb_gem_hit_vmoment, &b_Ndata_bb_gem_hit_vmoment);
   fChain->SetBranchAddress("bb.gem.hit.vmoment", bb_gem_hit_vmoment, &b_bb_gem_hit_vmoment);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vsigma", &Ndata_bb_gem_hit_vsigma, &b_Ndata_bb_gem_hit_vsigma);
   fChain->SetBranchAddress("bb.gem.hit.vsigma", bb_gem_hit_vsigma, &b_bb_gem_hit_vsigma);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstriphi", &Ndata_bb_gem_hit_vstriphi, &b_Ndata_bb_gem_hit_vstriphi);
   fChain->SetBranchAddress("bb.gem.hit.vstriphi", bb_gem_hit_vstriphi, &b_bb_gem_hit_vstriphi);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstriplo", &Ndata_bb_gem_hit_vstriplo, &b_Ndata_bb_gem_hit_vstriplo);
   fChain->SetBranchAddress("bb.gem.hit.vstriplo", bb_gem_hit_vstriplo, &b_bb_gem_hit_vstriplo);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.vstripmax", &Ndata_bb_gem_hit_vstripmax, &b_Ndata_bb_gem_hit_vstripmax);
   fChain->SetBranchAddress("bb.gem.hit.vstripmax", bb_gem_hit_vstripmax, &b_bb_gem_hit_vstripmax);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.xglobal", &Ndata_bb_gem_hit_xglobal, &b_Ndata_bb_gem_hit_xglobal);
   fChain->SetBranchAddress("bb.gem.hit.xglobal", bb_gem_hit_xglobal, &b_bb_gem_hit_xglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.xlocal", &Ndata_bb_gem_hit_xlocal, &b_Ndata_bb_gem_hit_xlocal);
   fChain->SetBranchAddress("bb.gem.hit.xlocal", bb_gem_hit_xlocal, &b_bb_gem_hit_xlocal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.yglobal", &Ndata_bb_gem_hit_yglobal, &b_Ndata_bb_gem_hit_yglobal);
   fChain->SetBranchAddress("bb.gem.hit.yglobal", bb_gem_hit_yglobal, &b_bb_gem_hit_yglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.ylocal", &Ndata_bb_gem_hit_ylocal, &b_Ndata_bb_gem_hit_ylocal);
   fChain->SetBranchAddress("bb.gem.hit.ylocal", bb_gem_hit_ylocal, &b_bb_gem_hit_ylocal);
   fChain->SetBranchAddress("Ndata.bb.gem.hit.zglobal", &Ndata_bb_gem_hit_zglobal, &b_Ndata_bb_gem_hit_zglobal);
   fChain->SetBranchAddress("bb.gem.hit.zglobal", bb_gem_hit_zglobal, &b_bb_gem_hit_zglobal);
   fChain->SetBranchAddress("Ndata.bb.gem.n2Dhit_layer", &Ndata_bb_gem_n2Dhit_layer, &b_Ndata_bb_gem_n2Dhit_layer);
   fChain->SetBranchAddress("bb.gem.n2Dhit_layer", bb_gem_n2Dhit_layer, &b_bb_gem_n2Dhit_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nclustu_layer", &Ndata_bb_gem_nclustu_layer, &b_Ndata_bb_gem_nclustu_layer);
   fChain->SetBranchAddress("bb.gem.nclustu_layer", bb_gem_nclustu_layer, &b_bb_gem_nclustu_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nclustv_layer", &Ndata_bb_gem_nclustv_layer, &b_Ndata_bb_gem_nclustv_layer);
   fChain->SetBranchAddress("bb.gem.nclustv_layer", bb_gem_nclustv_layer, &b_bb_gem_nclustv_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nstripsu_layer", &Ndata_bb_gem_nstripsu_layer, &b_Ndata_bb_gem_nstripsu_layer);
   fChain->SetBranchAddress("bb.gem.nstripsu_layer", bb_gem_nstripsu_layer, &b_bb_gem_nstripsu_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.nstripsv_layer", &Ndata_bb_gem_nstripsv_layer, &b_Ndata_bb_gem_nstripsv_layer);
   fChain->SetBranchAddress("bb.gem.nstripsv_layer", bb_gem_nstripsv_layer, &b_bb_gem_nstripsv_layer);
   fChain->SetBranchAddress("Ndata.bb.gem.track.chi2ndf", &Ndata_bb_gem_track_chi2ndf, &b_Ndata_bb_gem_track_chi2ndf);
   fChain->SetBranchAddress("bb.gem.track.chi2ndf", bb_gem_track_chi2ndf, &b_bb_gem_track_chi2ndf);
   fChain->SetBranchAddress("Ndata.bb.gem.track.nhits", &Ndata_bb_gem_track_nhits, &b_Ndata_bb_gem_track_nhits);
   fChain->SetBranchAddress("bb.gem.track.nhits", bb_gem_track_nhits, &b_bb_gem_track_nhits);
   fChain->SetBranchAddress("Ndata.bb.gem.track.x", &Ndata_bb_gem_track_x, &b_Ndata_bb_gem_track_x);
   fChain->SetBranchAddress("bb.gem.track.x", bb_gem_track_x, &b_bb_gem_track_x);
   fChain->SetBranchAddress("Ndata.bb.gem.track.xp", &Ndata_bb_gem_track_xp, &b_Ndata_bb_gem_track_xp);
   fChain->SetBranchAddress("bb.gem.track.xp", bb_gem_track_xp, &b_bb_gem_track_xp);
   fChain->SetBranchAddress("Ndata.bb.gem.track.y", &Ndata_bb_gem_track_y, &b_Ndata_bb_gem_track_y);
   fChain->SetBranchAddress("bb.gem.track.y", bb_gem_track_y, &b_bb_gem_track_y);
   fChain->SetBranchAddress("Ndata.bb.gem.track.yp", &Ndata_bb_gem_track_yp, &b_Ndata_bb_gem_track_yp);
   fChain->SetBranchAddress("bb.gem.track.yp", bb_gem_track_yp, &b_bb_gem_track_yp);
   fChain->SetBranchAddress("Ndata.bb.grinch.clus.adc", &Ndata_bb_grinch_clus_adc, &b_Ndata_bb_grinch_clus_adc);
   fChain->SetBranchAddress("bb.grinch.clus.adc", bb_grinch_clus_adc, &b_bb_grinch_clus_adc);
   fChain->SetBranchAddress("Ndata.bb.grinch.clus.size", &Ndata_bb_grinch_clus_size, &b_Ndata_bb_grinch_clus_size);
   fChain->SetBranchAddress("bb.grinch.clus.size", bb_grinch_clus_size, &b_bb_grinch_clus_size);
   fChain->SetBranchAddress("Ndata.bb.grinch.clus.tf_mean", &Ndata_bb_grinch_clus_tf_mean, &b_Ndata_bb_grinch_clus_tf_mean);
   fChain->SetBranchAddress("bb.grinch.clus.tf_mean", bb_grinch_clus_tf_mean, &b_bb_grinch_clus_tf_mean);
   fChain->SetBranchAddress("Ndata.bb.grinch.clus.tr_mean", &Ndata_bb_grinch_clus_tr_mean, &b_Ndata_bb_grinch_clus_tr_mean);
   fChain->SetBranchAddress("bb.grinch.clus.tr_mean", bb_grinch_clus_tr_mean, &b_bb_grinch_clus_tr_mean);
   fChain->SetBranchAddress("Ndata.bb.grinch.clus.x_mean", &Ndata_bb_grinch_clus_x_mean, &b_Ndata_bb_grinch_clus_x_mean);
   fChain->SetBranchAddress("bb.grinch.clus.x_mean", bb_grinch_clus_x_mean, &b_bb_grinch_clus_x_mean);
   fChain->SetBranchAddress("Ndata.bb.grinch.clus.y_mean", &Ndata_bb_grinch_clus_y_mean, &b_Ndata_bb_grinch_clus_y_mean);
   fChain->SetBranchAddress("bb.grinch.clus.y_mean", bb_grinch_clus_y_mean, &b_bb_grinch_clus_y_mean);
   fChain->SetBranchAddress("Ndata.bb.grinch.hit.adc", &Ndata_bb_grinch_hit_adc, &b_Ndata_bb_grinch_hit_adc);
   fChain->SetBranchAddress("bb.grinch.hit.adc", bb_grinch_hit_adc, &b_bb_grinch_hit_adc);
   fChain->SetBranchAddress("Ndata.bb.grinch.hit.col", &Ndata_bb_grinch_hit_col, &b_Ndata_bb_grinch_hit_col);
   fChain->SetBranchAddress("bb.grinch.hit.col", bb_grinch_hit_col, &b_bb_grinch_hit_col);
   fChain->SetBranchAddress("Ndata.bb.grinch.hit.pmtnum", &Ndata_bb_grinch_hit_pmtnum, &b_Ndata_bb_grinch_hit_pmtnum);
   fChain->SetBranchAddress("bb.grinch.hit.pmtnum", bb_grinch_hit_pmtnum, &b_bb_grinch_hit_pmtnum);
   fChain->SetBranchAddress("Ndata.bb.grinch.hit.row", &Ndata_bb_grinch_hit_row, &b_Ndata_bb_grinch_hit_row);
   fChain->SetBranchAddress("bb.grinch.hit.row", bb_grinch_hit_row, &b_bb_grinch_hit_row);
   fChain->SetBranchAddress("Ndata.bb.grinch.hit.tdc_f", &Ndata_bb_grinch_hit_tdc_f, &b_Ndata_bb_grinch_hit_tdc_f);
   fChain->SetBranchAddress("bb.grinch.hit.tdc_f", bb_grinch_hit_tdc_f, &b_bb_grinch_hit_tdc_f);
   fChain->SetBranchAddress("Ndata.bb.grinch.hit.tdc_r", &Ndata_bb_grinch_hit_tdc_r, &b_Ndata_bb_grinch_hit_tdc_r);
   fChain->SetBranchAddress("bb.grinch.hit.tdc_r", bb_grinch_hit_tdc_r, &b_bb_grinch_hit_tdc_r);
   fChain->SetBranchAddress("Ndata.bb.grinch.hit.xhit", &Ndata_bb_grinch_hit_xhit, &b_Ndata_bb_grinch_hit_xhit);
   fChain->SetBranchAddress("bb.grinch.hit.xhit", bb_grinch_hit_xhit, &b_bb_grinch_hit_xhit);
   fChain->SetBranchAddress("Ndata.bb.grinch.hit.yhit", &Ndata_bb_grinch_hit_yhit, &b_Ndata_bb_grinch_hit_yhit);
   fChain->SetBranchAddress("bb.grinch.hit.yhit", bb_grinch_hit_yhit, &b_bb_grinch_hit_yhit);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.L.le", &Ndata_bb_hodo_bar_tdc_L_le, &b_Ndata_bb_hodo_bar_tdc_L_le);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.L.le", bb_hodo_bar_tdc_L_le, &b_bb_hodo_bar_tdc_L_le);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.L.leW", &Ndata_bb_hodo_bar_tdc_L_leW, &b_Ndata_bb_hodo_bar_tdc_L_leW);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.L.leW", bb_hodo_bar_tdc_L_leW, &b_bb_hodo_bar_tdc_L_leW);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.L.te", &Ndata_bb_hodo_bar_tdc_L_te, &b_Ndata_bb_hodo_bar_tdc_L_te);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.L.te", bb_hodo_bar_tdc_L_te, &b_bb_hodo_bar_tdc_L_te);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.L.teW", &Ndata_bb_hodo_bar_tdc_L_teW, &b_Ndata_bb_hodo_bar_tdc_L_teW);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.L.teW", bb_hodo_bar_tdc_L_teW, &b_bb_hodo_bar_tdc_L_teW);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.L.tot", &Ndata_bb_hodo_bar_tdc_L_tot, &b_Ndata_bb_hodo_bar_tdc_L_tot);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.L.tot", bb_hodo_bar_tdc_L_tot, &b_bb_hodo_bar_tdc_L_tot);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.L.totW", &Ndata_bb_hodo_bar_tdc_L_totW, &b_Ndata_bb_hodo_bar_tdc_L_totW);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.L.totW", bb_hodo_bar_tdc_L_totW, &b_bb_hodo_bar_tdc_L_totW);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.R.le", &Ndata_bb_hodo_bar_tdc_R_le, &b_Ndata_bb_hodo_bar_tdc_R_le);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.R.le", bb_hodo_bar_tdc_R_le, &b_bb_hodo_bar_tdc_R_le);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.R.leW", &Ndata_bb_hodo_bar_tdc_R_leW, &b_Ndata_bb_hodo_bar_tdc_R_leW);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.R.leW", bb_hodo_bar_tdc_R_leW, &b_bb_hodo_bar_tdc_R_leW);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.R.te", &Ndata_bb_hodo_bar_tdc_R_te, &b_Ndata_bb_hodo_bar_tdc_R_te);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.R.te", bb_hodo_bar_tdc_R_te, &b_bb_hodo_bar_tdc_R_te);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.R.teW", &Ndata_bb_hodo_bar_tdc_R_teW, &b_Ndata_bb_hodo_bar_tdc_R_teW);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.R.teW", bb_hodo_bar_tdc_R_teW, &b_bb_hodo_bar_tdc_R_teW);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.R.tot", &Ndata_bb_hodo_bar_tdc_R_tot, &b_Ndata_bb_hodo_bar_tdc_R_tot);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.R.tot", bb_hodo_bar_tdc_R_tot, &b_bb_hodo_bar_tdc_R_tot);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.R.totW", &Ndata_bb_hodo_bar_tdc_R_totW, &b_Ndata_bb_hodo_bar_tdc_R_totW);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.R.totW", bb_hodo_bar_tdc_R_totW, &b_bb_hodo_bar_tdc_R_totW);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.id", &Ndata_bb_hodo_bar_tdc_id, &b_Ndata_bb_hodo_bar_tdc_id);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.id", bb_hodo_bar_tdc_id, &b_bb_hodo_bar_tdc_id);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.meantime", &Ndata_bb_hodo_bar_tdc_meantime, &b_Ndata_bb_hodo_bar_tdc_meantime);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.meantime", bb_hodo_bar_tdc_meantime, &b_bb_hodo_bar_tdc_meantime);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.timediff", &Ndata_bb_hodo_bar_tdc_timediff, &b_Ndata_bb_hodo_bar_tdc_timediff);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.timediff", bb_hodo_bar_tdc_timediff, &b_bb_hodo_bar_tdc_timediff);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.timehitpos", &Ndata_bb_hodo_bar_tdc_timehitpos, &b_Ndata_bb_hodo_bar_tdc_timehitpos);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.timehitpos", bb_hodo_bar_tdc_timehitpos, &b_bb_hodo_bar_tdc_timehitpos);
   fChain->SetBranchAddress("Ndata.bb.hodo.bar.tdc.vpos", &Ndata_bb_hodo_bar_tdc_vpos, &b_Ndata_bb_hodo_bar_tdc_vpos);
   fChain->SetBranchAddress("bb.hodo.bar.tdc.vpos", bb_hodo_bar_tdc_vpos, &b_bb_hodo_bar_tdc_vpos);
   fChain->SetBranchAddress("Ndata.bb.hodo.clus.size", &Ndata_bb_hodo_clus_size, &b_Ndata_bb_hodo_clus_size);
   fChain->SetBranchAddress("bb.hodo.clus.size", bb_hodo_clus_size, &b_bb_hodo_clus_size);
   fChain->SetBranchAddress("Ndata.bb.hodo.clus.tmean", &Ndata_bb_hodo_clus_tmean, &b_Ndata_bb_hodo_clus_tmean);
   fChain->SetBranchAddress("bb.hodo.clus.tmean", bb_hodo_clus_tmean, &b_bb_hodo_clus_tmean);
   fChain->SetBranchAddress("Ndata.bb.hodo.clus.totmean", &Ndata_bb_hodo_clus_totmean, &b_Ndata_bb_hodo_clus_totmean);
   fChain->SetBranchAddress("bb.hodo.clus.totmean", bb_hodo_clus_totmean, &b_bb_hodo_clus_totmean);
   fChain->SetBranchAddress("Ndata.bb.hodo.clus.xmean", &Ndata_bb_hodo_clus_xmean, &b_Ndata_bb_hodo_clus_xmean);
   fChain->SetBranchAddress("bb.hodo.clus.xmean", bb_hodo_clus_xmean, &b_bb_hodo_clus_xmean);
   fChain->SetBranchAddress("Ndata.bb.hodo.clus.ymean", &Ndata_bb_hodo_clus_ymean, &b_Ndata_bb_hodo_clus_ymean);
   fChain->SetBranchAddress("bb.hodo.clus.ymean", bb_hodo_clus_ymean, &b_bb_hodo_clus_ymean);
   fChain->SetBranchAddress("Ndata.bb.hodo.tdc", &Ndata_bb_hodo_tdc, &b_Ndata_bb_hodo_tdc);
   fChain->SetBranchAddress("bb.hodo.tdc", bb_hodo_tdc, &b_bb_hodo_tdc);
   fChain->SetBranchAddress("Ndata.bb.hodo.tdc_mult", &Ndata_bb_hodo_tdc_mult, &b_Ndata_bb_hodo_tdc_mult);
   fChain->SetBranchAddress("bb.hodo.tdc_mult", bb_hodo_tdc_mult, &b_bb_hodo_tdc_mult);
   fChain->SetBranchAddress("Ndata.bb.hodo.tdc_te", &Ndata_bb_hodo_tdc_te, &b_Ndata_bb_hodo_tdc_te);
   fChain->SetBranchAddress("bb.hodo.tdc_te", bb_hodo_tdc_te, &b_bb_hodo_tdc_te);
   fChain->SetBranchAddress("Ndata.bb.hodo.tdc_tot", &Ndata_bb_hodo_tdc_tot, &b_Ndata_bb_hodo_tdc_tot);
   fChain->SetBranchAddress("bb.hodo.tdc_tot", bb_hodo_tdc_tot, &b_bb_hodo_tdc_tot);
   fChain->SetBranchAddress("Ndata.bb.hodo.tdccol", &Ndata_bb_hodo_tdccol, &b_Ndata_bb_hodo_tdccol);
   fChain->SetBranchAddress("bb.hodo.tdccol", bb_hodo_tdccol, &b_bb_hodo_tdccol);
   fChain->SetBranchAddress("Ndata.bb.hodo.tdcelemID", &Ndata_bb_hodo_tdcelemID, &b_Ndata_bb_hodo_tdcelemID);
   fChain->SetBranchAddress("bb.hodo.tdcelemID", bb_hodo_tdcelemID, &b_bb_hodo_tdcelemID);
   fChain->SetBranchAddress("Ndata.bb.hodo.tdclayer", &Ndata_bb_hodo_tdclayer, &b_Ndata_bb_hodo_tdclayer);
   fChain->SetBranchAddress("bb.hodo.tdclayer", bb_hodo_tdclayer, &b_bb_hodo_tdclayer);
   fChain->SetBranchAddress("Ndata.bb.hodo.tdcrow", &Ndata_bb_hodo_tdcrow, &b_Ndata_bb_hodo_tdcrow);
   fChain->SetBranchAddress("bb.hodo.tdcrow", bb_hodo_tdcrow, &b_bb_hodo_tdcrow);
   fChain->SetBranchAddress("Ndata.bb.prob_e", &Ndata_bb_prob_e, &b_Ndata_bb_prob_e);
   fChain->SetBranchAddress("bb.prob_e", &bb_prob_e, &b_bb_prob_e);
   fChain->SetBranchAddress("Ndata.bb.prob_pi", &Ndata_bb_prob_pi, &b_Ndata_bb_prob_pi);
   fChain->SetBranchAddress("bb.prob_pi", &bb_prob_pi, &b_bb_prob_pi);
   fChain->SetBranchAddress("Ndata.bb.ps.a", &Ndata_bb_ps_a, &b_Ndata_bb_ps_a);
   fChain->SetBranchAddress("bb.ps.a", bb_ps_a, &b_bb_ps_a);
   fChain->SetBranchAddress("Ndata.bb.ps.a_amp", &Ndata_bb_ps_a_amp, &b_Ndata_bb_ps_a_amp);
   fChain->SetBranchAddress("bb.ps.a_amp", bb_ps_a_amp, &b_bb_ps_a_amp);
   fChain->SetBranchAddress("Ndata.bb.ps.a_amp_p", &Ndata_bb_ps_a_amp_p, &b_Ndata_bb_ps_a_amp_p);
   fChain->SetBranchAddress("bb.ps.a_amp_p", bb_ps_a_amp_p, &b_bb_ps_a_amp_p);
   fChain->SetBranchAddress("Ndata.bb.ps.a_c", &Ndata_bb_ps_a_c, &b_Ndata_bb_ps_a_c);
   fChain->SetBranchAddress("bb.ps.a_c", bb_ps_a_c, &b_bb_ps_a_c);
   fChain->SetBranchAddress("Ndata.bb.ps.a_mult", &Ndata_bb_ps_a_mult, &b_Ndata_bb_ps_a_mult);
   fChain->SetBranchAddress("bb.ps.a_mult", &bb_ps_a_mult, &b_bb_ps_a_mult);
   fChain->SetBranchAddress("Ndata.bb.ps.a_p", &Ndata_bb_ps_a_p, &b_Ndata_bb_ps_a_p);
   fChain->SetBranchAddress("bb.ps.a_p", bb_ps_a_p, &b_bb_ps_a_p);
   fChain->SetBranchAddress("Ndata.bb.ps.a_time", &Ndata_bb_ps_a_time, &b_Ndata_bb_ps_a_time);
   fChain->SetBranchAddress("bb.ps.a_time", bb_ps_a_time, &b_bb_ps_a_time);
   fChain->SetBranchAddress("Ndata.bb.ps.adccol", &Ndata_bb_ps_adccol, &b_Ndata_bb_ps_adccol);
   fChain->SetBranchAddress("bb.ps.adccol", bb_ps_adccol, &b_bb_ps_adccol);
   fChain->SetBranchAddress("Ndata.bb.ps.adcelemID", &Ndata_bb_ps_adcelemID, &b_Ndata_bb_ps_adcelemID);
   fChain->SetBranchAddress("bb.ps.adcelemID", bb_ps_adcelemID, &b_bb_ps_adcelemID);
   fChain->SetBranchAddress("Ndata.bb.ps.adclayer", &Ndata_bb_ps_adclayer, &b_Ndata_bb_ps_adclayer);
   fChain->SetBranchAddress("bb.ps.adclayer", bb_ps_adclayer, &b_bb_ps_adclayer);
   fChain->SetBranchAddress("Ndata.bb.ps.adcrow", &Ndata_bb_ps_adcrow, &b_Ndata_bb_ps_adcrow);
   fChain->SetBranchAddress("bb.ps.adcrow", bb_ps_adcrow, &b_bb_ps_adcrow);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.col", &Ndata_bb_ps_clus_col, &b_Ndata_bb_ps_clus_col);
   fChain->SetBranchAddress("bb.ps.clus.col", bb_ps_clus_col, &b_bb_ps_clus_col);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.e", &Ndata_bb_ps_clus_e, &b_Ndata_bb_ps_clus_e);
   fChain->SetBranchAddress("bb.ps.clus.e", bb_ps_clus_e, &b_bb_ps_clus_e);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.e_c", &Ndata_bb_ps_clus_e_c, &b_Ndata_bb_ps_clus_e_c);
   fChain->SetBranchAddress("bb.ps.clus.e_c", bb_ps_clus_e_c, &b_bb_ps_clus_e_c);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.eblk", &Ndata_bb_ps_clus_eblk, &b_Ndata_bb_ps_clus_eblk);
   fChain->SetBranchAddress("bb.ps.clus.eblk", bb_ps_clus_eblk, &b_bb_ps_clus_eblk);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.eblk_c", &Ndata_bb_ps_clus_eblk_c, &b_Ndata_bb_ps_clus_eblk_c);
   fChain->SetBranchAddress("bb.ps.clus.eblk_c", bb_ps_clus_eblk_c, &b_bb_ps_clus_eblk_c);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.id", &Ndata_bb_ps_clus_id, &b_Ndata_bb_ps_clus_id);
   fChain->SetBranchAddress("bb.ps.clus.id", bb_ps_clus_id, &b_bb_ps_clus_id);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.nblk", &Ndata_bb_ps_clus_nblk, &b_Ndata_bb_ps_clus_nblk);
   fChain->SetBranchAddress("bb.ps.clus.nblk", bb_ps_clus_nblk, &b_bb_ps_clus_nblk);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.row", &Ndata_bb_ps_clus_row, &b_Ndata_bb_ps_clus_row);
   fChain->SetBranchAddress("bb.ps.clus.row", bb_ps_clus_row, &b_bb_ps_clus_row);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.x", &Ndata_bb_ps_clus_x, &b_Ndata_bb_ps_clus_x);
   fChain->SetBranchAddress("bb.ps.clus.x", bb_ps_clus_x, &b_bb_ps_clus_x);
   fChain->SetBranchAddress("Ndata.bb.ps.clus.y", &Ndata_bb_ps_clus_y, &b_Ndata_bb_ps_clus_y);
   fChain->SetBranchAddress("bb.ps.clus.y", bb_ps_clus_y, &b_bb_ps_clus_y);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.col", &Ndata_bb_ps_clus_blk_col, &b_Ndata_bb_ps_clus_blk_col);
   fChain->SetBranchAddress("bb.ps.clus_blk.col", bb_ps_clus_blk_col, &b_bb_ps_clus_blk_col);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.e", &Ndata_bb_ps_clus_blk_e, &b_Ndata_bb_ps_clus_blk_e);
   fChain->SetBranchAddress("bb.ps.clus_blk.e", bb_ps_clus_blk_e, &b_bb_ps_clus_blk_e);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.e_c", &Ndata_bb_ps_clus_blk_e_c, &b_Ndata_bb_ps_clus_blk_e_c);
   fChain->SetBranchAddress("bb.ps.clus_blk.e_c", &bb_ps_clus_blk_e_c, &b_bb_ps_clus_blk_e_c);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.id", &Ndata_bb_ps_clus_blk_id, &b_Ndata_bb_ps_clus_blk_id);
   fChain->SetBranchAddress("bb.ps.clus_blk.id", bb_ps_clus_blk_id, &b_bb_ps_clus_blk_id);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.row", &Ndata_bb_ps_clus_blk_row, &b_Ndata_bb_ps_clus_blk_row);
   fChain->SetBranchAddress("bb.ps.clus_blk.row", bb_ps_clus_blk_row, &b_bb_ps_clus_blk_row);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.x", &Ndata_bb_ps_clus_blk_x, &b_Ndata_bb_ps_clus_blk_x);
   fChain->SetBranchAddress("bb.ps.clus_blk.x", bb_ps_clus_blk_x, &b_bb_ps_clus_blk_x);
   fChain->SetBranchAddress("Ndata.bb.ps.clus_blk.y", &Ndata_bb_ps_clus_blk_y, &b_Ndata_bb_ps_clus_blk_y);
   fChain->SetBranchAddress("bb.ps.clus_blk.y", bb_ps_clus_blk_y, &b_bb_ps_clus_blk_y);
   fChain->SetBranchAddress("Ndata.bb.ps.e_res", &Ndata_bb_ps_e_res, &b_Ndata_bb_ps_e_res);
   fChain->SetBranchAddress("bb.ps.e_res", bb_ps_e_res, &b_bb_ps_e_res);
   fChain->SetBranchAddress("Ndata.bb.ps.goodblock.atime", &Ndata_bb_ps_goodblock_atime, &b_Ndata_bb_ps_goodblock_atime);
   fChain->SetBranchAddress("bb.ps.goodblock.atime", bb_ps_goodblock_atime, &b_bb_ps_goodblock_atime);
   fChain->SetBranchAddress("Ndata.bb.ps.goodblock.col", &Ndata_bb_ps_goodblock_col, &b_Ndata_bb_ps_goodblock_col);
   fChain->SetBranchAddress("bb.ps.goodblock.col", bb_ps_goodblock_col, &b_bb_ps_goodblock_col);
   fChain->SetBranchAddress("Ndata.bb.ps.goodblock.e", &Ndata_bb_ps_goodblock_e, &b_Ndata_bb_ps_goodblock_e);
   fChain->SetBranchAddress("bb.ps.goodblock.e", bb_ps_goodblock_e, &b_bb_ps_goodblock_e);
   fChain->SetBranchAddress("Ndata.bb.ps.goodblock.id", &Ndata_bb_ps_goodblock_id, &b_Ndata_bb_ps_goodblock_id);
   fChain->SetBranchAddress("bb.ps.goodblock.id", bb_ps_goodblock_id, &b_bb_ps_goodblock_id);
   fChain->SetBranchAddress("Ndata.bb.ps.goodblock.row", &Ndata_bb_ps_goodblock_row, &b_Ndata_bb_ps_goodblock_row);
   fChain->SetBranchAddress("bb.ps.goodblock.row", bb_ps_goodblock_row, &b_bb_ps_goodblock_row);
   fChain->SetBranchAddress("Ndata.bb.ps.goodblock.x", &Ndata_bb_ps_goodblock_x, &b_Ndata_bb_ps_goodblock_x);
   fChain->SetBranchAddress("bb.ps.goodblock.x", bb_ps_goodblock_x, &b_bb_ps_goodblock_x);
   fChain->SetBranchAddress("Ndata.bb.ps.goodblock.y", &Ndata_bb_ps_goodblock_y, &b_Ndata_bb_ps_goodblock_y);
   fChain->SetBranchAddress("bb.ps.goodblock.y", bb_ps_goodblock_y, &b_bb_ps_goodblock_y);
   fChain->SetBranchAddress("Ndata.bb.ps.ped", &Ndata_bb_ps_ped, &b_Ndata_bb_ps_ped);
   fChain->SetBranchAddress("bb.ps.ped", bb_ps_ped, &b_bb_ps_ped);
   fChain->SetBranchAddress("Ndata.bb.ps.x_res", &Ndata_bb_ps_x_res, &b_Ndata_bb_ps_x_res);
   fChain->SetBranchAddress("bb.ps.x_res", bb_ps_x_res, &b_bb_ps_x_res);
   fChain->SetBranchAddress("Ndata.bb.ps.y_res", &Ndata_bb_ps_y_res, &b_Ndata_bb_ps_y_res);
   fChain->SetBranchAddress("bb.ps.y_res", bb_ps_y_res, &b_bb_ps_y_res);
   fChain->SetBranchAddress("Ndata.bb.sh.a", &Ndata_bb_sh_a, &b_Ndata_bb_sh_a);
   fChain->SetBranchAddress("bb.sh.a", bb_sh_a, &b_bb_sh_a);
   fChain->SetBranchAddress("Ndata.bb.sh.a_amp", &Ndata_bb_sh_a_amp, &b_Ndata_bb_sh_a_amp);
   fChain->SetBranchAddress("bb.sh.a_amp", bb_sh_a_amp, &b_bb_sh_a_amp);
   fChain->SetBranchAddress("Ndata.bb.sh.a_amp_p", &Ndata_bb_sh_a_amp_p, &b_Ndata_bb_sh_a_amp_p);
   fChain->SetBranchAddress("bb.sh.a_amp_p", bb_sh_a_amp_p, &b_bb_sh_a_amp_p);
   fChain->SetBranchAddress("Ndata.bb.sh.a_c", &Ndata_bb_sh_a_c, &b_Ndata_bb_sh_a_c);
   fChain->SetBranchAddress("bb.sh.a_c", bb_sh_a_c, &b_bb_sh_a_c);
   fChain->SetBranchAddress("Ndata.bb.sh.a_mult", &Ndata_bb_sh_a_mult, &b_Ndata_bb_sh_a_mult);
   fChain->SetBranchAddress("bb.sh.a_mult", &bb_sh_a_mult, &b_bb_sh_a_mult);
   fChain->SetBranchAddress("Ndata.bb.sh.a_p", &Ndata_bb_sh_a_p, &b_Ndata_bb_sh_a_p);
   fChain->SetBranchAddress("bb.sh.a_p", bb_sh_a_p, &b_bb_sh_a_p);
   fChain->SetBranchAddress("Ndata.bb.sh.a_time", &Ndata_bb_sh_a_time, &b_Ndata_bb_sh_a_time);
   fChain->SetBranchAddress("bb.sh.a_time", bb_sh_a_time, &b_bb_sh_a_time);
   fChain->SetBranchAddress("Ndata.bb.sh.adccol", &Ndata_bb_sh_adccol, &b_Ndata_bb_sh_adccol);
   fChain->SetBranchAddress("bb.sh.adccol", bb_sh_adccol, &b_bb_sh_adccol);
   fChain->SetBranchAddress("Ndata.bb.sh.adcelemID", &Ndata_bb_sh_adcelemID, &b_Ndata_bb_sh_adcelemID);
   fChain->SetBranchAddress("bb.sh.adcelemID", bb_sh_adcelemID, &b_bb_sh_adcelemID);
   fChain->SetBranchAddress("Ndata.bb.sh.adclayer", &Ndata_bb_sh_adclayer, &b_Ndata_bb_sh_adclayer);
   fChain->SetBranchAddress("bb.sh.adclayer", bb_sh_adclayer, &b_bb_sh_adclayer);
   fChain->SetBranchAddress("Ndata.bb.sh.adcrow", &Ndata_bb_sh_adcrow, &b_Ndata_bb_sh_adcrow);
   fChain->SetBranchAddress("bb.sh.adcrow", bb_sh_adcrow, &b_bb_sh_adcrow);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.col", &Ndata_bb_sh_clus_col, &b_Ndata_bb_sh_clus_col);
   fChain->SetBranchAddress("bb.sh.clus.col", bb_sh_clus_col, &b_bb_sh_clus_col);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.e", &Ndata_bb_sh_clus_e, &b_Ndata_bb_sh_clus_e);
   fChain->SetBranchAddress("bb.sh.clus.e", bb_sh_clus_e, &b_bb_sh_clus_e);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.e_c", &Ndata_bb_sh_clus_e_c, &b_Ndata_bb_sh_clus_e_c);
   fChain->SetBranchAddress("bb.sh.clus.e_c", bb_sh_clus_e_c, &b_bb_sh_clus_e_c);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.eblk", &Ndata_bb_sh_clus_eblk, &b_Ndata_bb_sh_clus_eblk);
   fChain->SetBranchAddress("bb.sh.clus.eblk", bb_sh_clus_eblk, &b_bb_sh_clus_eblk);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.eblk_c", &Ndata_bb_sh_clus_eblk_c, &b_Ndata_bb_sh_clus_eblk_c);
   fChain->SetBranchAddress("bb.sh.clus.eblk_c", bb_sh_clus_eblk_c, &b_bb_sh_clus_eblk_c);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.id", &Ndata_bb_sh_clus_id, &b_Ndata_bb_sh_clus_id);
   fChain->SetBranchAddress("bb.sh.clus.id", bb_sh_clus_id, &b_bb_sh_clus_id);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.nblk", &Ndata_bb_sh_clus_nblk, &b_Ndata_bb_sh_clus_nblk);
   fChain->SetBranchAddress("bb.sh.clus.nblk", bb_sh_clus_nblk, &b_bb_sh_clus_nblk);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.row", &Ndata_bb_sh_clus_row, &b_Ndata_bb_sh_clus_row);
   fChain->SetBranchAddress("bb.sh.clus.row", bb_sh_clus_row, &b_bb_sh_clus_row);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.x", &Ndata_bb_sh_clus_x, &b_Ndata_bb_sh_clus_x);
   fChain->SetBranchAddress("bb.sh.clus.x", bb_sh_clus_x, &b_bb_sh_clus_x);
   fChain->SetBranchAddress("Ndata.bb.sh.clus.y", &Ndata_bb_sh_clus_y, &b_Ndata_bb_sh_clus_y);
   fChain->SetBranchAddress("bb.sh.clus.y", bb_sh_clus_y, &b_bb_sh_clus_y);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.col", &Ndata_bb_sh_clus_blk_col, &b_Ndata_bb_sh_clus_blk_col);
   fChain->SetBranchAddress("bb.sh.clus_blk.col", bb_sh_clus_blk_col, &b_bb_sh_clus_blk_col);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.e", &Ndata_bb_sh_clus_blk_e, &b_Ndata_bb_sh_clus_blk_e);
   fChain->SetBranchAddress("bb.sh.clus_blk.e", bb_sh_clus_blk_e, &b_bb_sh_clus_blk_e);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.e_c", &Ndata_bb_sh_clus_blk_e_c, &b_Ndata_bb_sh_clus_blk_e_c);
   fChain->SetBranchAddress("bb.sh.clus_blk.e_c", &bb_sh_clus_blk_e_c, &b_bb_sh_clus_blk_e_c);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.id", &Ndata_bb_sh_clus_blk_id, &b_Ndata_bb_sh_clus_blk_id);
   fChain->SetBranchAddress("bb.sh.clus_blk.id", bb_sh_clus_blk_id, &b_bb_sh_clus_blk_id);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.row", &Ndata_bb_sh_clus_blk_row, &b_Ndata_bb_sh_clus_blk_row);
   fChain->SetBranchAddress("bb.sh.clus_blk.row", bb_sh_clus_blk_row, &b_bb_sh_clus_blk_row);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.x", &Ndata_bb_sh_clus_blk_x, &b_Ndata_bb_sh_clus_blk_x);
   fChain->SetBranchAddress("bb.sh.clus_blk.x", bb_sh_clus_blk_x, &b_bb_sh_clus_blk_x);
   fChain->SetBranchAddress("Ndata.bb.sh.clus_blk.y", &Ndata_bb_sh_clus_blk_y, &b_Ndata_bb_sh_clus_blk_y);
   fChain->SetBranchAddress("bb.sh.clus_blk.y", bb_sh_clus_blk_y, &b_bb_sh_clus_blk_y);
   fChain->SetBranchAddress("Ndata.bb.sh.e_res", &Ndata_bb_sh_e_res, &b_Ndata_bb_sh_e_res);
   fChain->SetBranchAddress("bb.sh.e_res", bb_sh_e_res, &b_bb_sh_e_res);
   fChain->SetBranchAddress("Ndata.bb.sh.goodblock.atime", &Ndata_bb_sh_goodblock_atime, &b_Ndata_bb_sh_goodblock_atime);
   fChain->SetBranchAddress("bb.sh.goodblock.atime", bb_sh_goodblock_atime, &b_bb_sh_goodblock_atime);
   fChain->SetBranchAddress("Ndata.bb.sh.goodblock.col", &Ndata_bb_sh_goodblock_col, &b_Ndata_bb_sh_goodblock_col);
   fChain->SetBranchAddress("bb.sh.goodblock.col", bb_sh_goodblock_col, &b_bb_sh_goodblock_col);
   fChain->SetBranchAddress("Ndata.bb.sh.goodblock.e", &Ndata_bb_sh_goodblock_e, &b_Ndata_bb_sh_goodblock_e);
   fChain->SetBranchAddress("bb.sh.goodblock.e", bb_sh_goodblock_e, &b_bb_sh_goodblock_e);
   fChain->SetBranchAddress("Ndata.bb.sh.goodblock.id", &Ndata_bb_sh_goodblock_id, &b_Ndata_bb_sh_goodblock_id);
   fChain->SetBranchAddress("bb.sh.goodblock.id", bb_sh_goodblock_id, &b_bb_sh_goodblock_id);
   fChain->SetBranchAddress("Ndata.bb.sh.goodblock.row", &Ndata_bb_sh_goodblock_row, &b_Ndata_bb_sh_goodblock_row);
   fChain->SetBranchAddress("bb.sh.goodblock.row", bb_sh_goodblock_row, &b_bb_sh_goodblock_row);
   fChain->SetBranchAddress("Ndata.bb.sh.goodblock.x", &Ndata_bb_sh_goodblock_x, &b_Ndata_bb_sh_goodblock_x);
   fChain->SetBranchAddress("bb.sh.goodblock.x", bb_sh_goodblock_x, &b_bb_sh_goodblock_x);
   fChain->SetBranchAddress("Ndata.bb.sh.goodblock.y", &Ndata_bb_sh_goodblock_y, &b_Ndata_bb_sh_goodblock_y);
   fChain->SetBranchAddress("bb.sh.goodblock.y", bb_sh_goodblock_y, &b_bb_sh_goodblock_y);
   fChain->SetBranchAddress("Ndata.bb.sh.ped", &Ndata_bb_sh_ped, &b_Ndata_bb_sh_ped);
   fChain->SetBranchAddress("bb.sh.ped", bb_sh_ped, &b_bb_sh_ped);
   fChain->SetBranchAddress("Ndata.bb.sh.x_res", &Ndata_bb_sh_x_res, &b_Ndata_bb_sh_x_res);
   fChain->SetBranchAddress("bb.sh.x_res", bb_sh_x_res, &b_bb_sh_x_res);
   fChain->SetBranchAddress("Ndata.bb.sh.y_res", &Ndata_bb_sh_y_res, &b_Ndata_bb_sh_y_res);
   fChain->SetBranchAddress("bb.sh.y_res", bb_sh_y_res, &b_bb_sh_y_res);
   fChain->SetBranchAddress("Ndata.bb.tr.beta", &Ndata_bb_tr_beta, &b_Ndata_bb_tr_beta);
   fChain->SetBranchAddress("bb.tr.beta", bb_tr_beta, &b_bb_tr_beta);
   fChain->SetBranchAddress("Ndata.bb.tr.chi2", &Ndata_bb_tr_chi2, &b_Ndata_bb_tr_chi2);
   fChain->SetBranchAddress("bb.tr.chi2", bb_tr_chi2, &b_bb_tr_chi2);
   fChain->SetBranchAddress("Ndata.bb.tr.d_ph", &Ndata_bb_tr_d_ph, &b_Ndata_bb_tr_d_ph);
   fChain->SetBranchAddress("bb.tr.d_ph", bb_tr_d_ph, &b_bb_tr_d_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.d_th", &Ndata_bb_tr_d_th, &b_Ndata_bb_tr_d_th);
   fChain->SetBranchAddress("bb.tr.d_th", bb_tr_d_th, &b_bb_tr_d_th);
   fChain->SetBranchAddress("Ndata.bb.tr.d_x", &Ndata_bb_tr_d_x, &b_Ndata_bb_tr_d_x);
   fChain->SetBranchAddress("bb.tr.d_x", bb_tr_d_x, &b_bb_tr_d_x);
   fChain->SetBranchAddress("Ndata.bb.tr.d_y", &Ndata_bb_tr_d_y, &b_Ndata_bb_tr_d_y);
   fChain->SetBranchAddress("bb.tr.d_y", bb_tr_d_y, &b_bb_tr_d_y);
   fChain->SetBranchAddress("Ndata.bb.tr.dbeta", &Ndata_bb_tr_dbeta, &b_Ndata_bb_tr_dbeta);
   fChain->SetBranchAddress("bb.tr.dbeta", bb_tr_dbeta, &b_bb_tr_dbeta);
   fChain->SetBranchAddress("Ndata.bb.tr.dtime", &Ndata_bb_tr_dtime, &b_Ndata_bb_tr_dtime);
   fChain->SetBranchAddress("bb.tr.dtime", bb_tr_dtime, &b_bb_tr_dtime);
   fChain->SetBranchAddress("Ndata.bb.tr.flag", &Ndata_bb_tr_flag, &b_Ndata_bb_tr_flag);
   fChain->SetBranchAddress("bb.tr.flag", bb_tr_flag, &b_bb_tr_flag);
   fChain->SetBranchAddress("Ndata.bb.tr.ndof", &Ndata_bb_tr_ndof, &b_Ndata_bb_tr_ndof);
   fChain->SetBranchAddress("bb.tr.ndof", bb_tr_ndof, &b_bb_tr_ndof);
   fChain->SetBranchAddress("Ndata.bb.tr.p", &Ndata_bb_tr_p, &b_Ndata_bb_tr_p);
   fChain->SetBranchAddress("bb.tr.p", bb_tr_p, &b_bb_tr_p);
   fChain->SetBranchAddress("Ndata.bb.tr.pathl", &Ndata_bb_tr_pathl, &b_Ndata_bb_tr_pathl);
   fChain->SetBranchAddress("bb.tr.pathl", bb_tr_pathl, &b_bb_tr_pathl);
   fChain->SetBranchAddress("Ndata.bb.tr.ph", &Ndata_bb_tr_ph, &b_Ndata_bb_tr_ph);
   fChain->SetBranchAddress("bb.tr.ph", bb_tr_ph, &b_bb_tr_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.px", &Ndata_bb_tr_px, &b_Ndata_bb_tr_px);
   fChain->SetBranchAddress("bb.tr.px", bb_tr_px, &b_bb_tr_px);
   fChain->SetBranchAddress("Ndata.bb.tr.py", &Ndata_bb_tr_py, &b_Ndata_bb_tr_py);
   fChain->SetBranchAddress("bb.tr.py", bb_tr_py, &b_bb_tr_py);
   fChain->SetBranchAddress("Ndata.bb.tr.pz", &Ndata_bb_tr_pz, &b_Ndata_bb_tr_pz);
   fChain->SetBranchAddress("bb.tr.pz", bb_tr_pz, &b_bb_tr_pz);
   fChain->SetBranchAddress("Ndata.bb.tr.r_ph", &Ndata_bb_tr_r_ph, &b_Ndata_bb_tr_r_ph);
   fChain->SetBranchAddress("bb.tr.r_ph", bb_tr_r_ph, &b_bb_tr_r_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.r_th", &Ndata_bb_tr_r_th, &b_Ndata_bb_tr_r_th);
   fChain->SetBranchAddress("bb.tr.r_th", bb_tr_r_th, &b_bb_tr_r_th);
   fChain->SetBranchAddress("Ndata.bb.tr.r_x", &Ndata_bb_tr_r_x, &b_Ndata_bb_tr_r_x);
   fChain->SetBranchAddress("bb.tr.r_x", bb_tr_r_x, &b_bb_tr_r_x);
   fChain->SetBranchAddress("Ndata.bb.tr.r_y", &Ndata_bb_tr_r_y, &b_Ndata_bb_tr_r_y);
   fChain->SetBranchAddress("bb.tr.r_y", bb_tr_r_y, &b_bb_tr_r_y);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_dp", &Ndata_bb_tr_tg_dp, &b_Ndata_bb_tr_tg_dp);
   fChain->SetBranchAddress("bb.tr.tg_dp", bb_tr_tg_dp, &b_bb_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_ph", &Ndata_bb_tr_tg_ph, &b_Ndata_bb_tr_tg_ph);
   fChain->SetBranchAddress("bb.tr.tg_ph", bb_tr_tg_ph, &b_bb_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_th", &Ndata_bb_tr_tg_th, &b_Ndata_bb_tr_tg_th);
   fChain->SetBranchAddress("bb.tr.tg_th", bb_tr_tg_th, &b_bb_tr_tg_th);
   fChain->SetBranchAddress("Ndata.bb.tr.tg_y", &Ndata_bb_tr_tg_y, &b_Ndata_bb_tr_tg_y);
   fChain->SetBranchAddress("bb.tr.tg_y", bb_tr_tg_y, &b_bb_tr_tg_y);
   fChain->SetBranchAddress("Ndata.bb.tr.th", &Ndata_bb_tr_th, &b_Ndata_bb_tr_th);
   fChain->SetBranchAddress("bb.tr.th", bb_tr_th, &b_bb_tr_th);
   fChain->SetBranchAddress("Ndata.bb.tr.time", &Ndata_bb_tr_time, &b_Ndata_bb_tr_time);
   fChain->SetBranchAddress("bb.tr.time", bb_tr_time, &b_bb_tr_time);
   fChain->SetBranchAddress("Ndata.bb.tr.vx", &Ndata_bb_tr_vx, &b_Ndata_bb_tr_vx);
   fChain->SetBranchAddress("bb.tr.vx", bb_tr_vx, &b_bb_tr_vx);
   fChain->SetBranchAddress("Ndata.bb.tr.vy", &Ndata_bb_tr_vy, &b_Ndata_bb_tr_vy);
   fChain->SetBranchAddress("bb.tr.vy", bb_tr_vy, &b_bb_tr_vy);
   fChain->SetBranchAddress("Ndata.bb.tr.vz", &Ndata_bb_tr_vz, &b_Ndata_bb_tr_vz);
   fChain->SetBranchAddress("bb.tr.vz", bb_tr_vz, &b_bb_tr_vz);
   fChain->SetBranchAddress("Ndata.bb.tr.x", &Ndata_bb_tr_x, &b_Ndata_bb_tr_x);
   fChain->SetBranchAddress("bb.tr.x", bb_tr_x, &b_bb_tr_x);
   fChain->SetBranchAddress("Ndata.bb.tr.y", &Ndata_bb_tr_y, &b_Ndata_bb_tr_y);
   fChain->SetBranchAddress("bb.tr.y", bb_tr_y, &b_bb_tr_y);
   fChain->SetBranchAddress("Ndata.bb.x_bcp", &Ndata_bb_x_bcp, &b_Ndata_bb_x_bcp);
   fChain->SetBranchAddress("bb.x_bcp", bb_x_bcp, &b_bb_x_bcp);
   fChain->SetBranchAddress("Ndata.bb.x_fcp", &Ndata_bb_x_fcp, &b_Ndata_bb_x_fcp);
   fChain->SetBranchAddress("bb.x_fcp", bb_x_fcp, &b_bb_x_fcp);
   fChain->SetBranchAddress("Ndata.bb.y_bcp", &Ndata_bb_y_bcp, &b_Ndata_bb_y_bcp);
   fChain->SetBranchAddress("bb.y_bcp", bb_y_bcp, &b_bb_y_bcp);
   fChain->SetBranchAddress("Ndata.bb.y_fcp", &Ndata_bb_y_fcp, &b_Ndata_bb_y_fcp);
   fChain->SetBranchAddress("bb.y_fcp", bb_y_fcp, &b_bb_y_fcp);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a", &Ndata_sbs_hcal_a, &b_Ndata_sbs_hcal_a);
   fChain->SetBranchAddress("sbs.hcal.a", sbs_hcal_a, &b_sbs_hcal_a);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_amp", &Ndata_sbs_hcal_a_amp, &b_Ndata_sbs_hcal_a_amp);
   fChain->SetBranchAddress("sbs.hcal.a_amp", sbs_hcal_a_amp, &b_sbs_hcal_a_amp);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_amp_p", &Ndata_sbs_hcal_a_amp_p, &b_Ndata_sbs_hcal_a_amp_p);
   fChain->SetBranchAddress("sbs.hcal.a_amp_p", sbs_hcal_a_amp_p, &b_sbs_hcal_a_amp_p);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_c", &Ndata_sbs_hcal_a_c, &b_Ndata_sbs_hcal_a_c);
   fChain->SetBranchAddress("sbs.hcal.a_c", sbs_hcal_a_c, &b_sbs_hcal_a_c);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_mult", &Ndata_sbs_hcal_a_mult, &b_Ndata_sbs_hcal_a_mult);
   fChain->SetBranchAddress("sbs.hcal.a_mult", &sbs_hcal_a_mult, &b_sbs_hcal_a_mult);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_p", &Ndata_sbs_hcal_a_p, &b_Ndata_sbs_hcal_a_p);
   fChain->SetBranchAddress("sbs.hcal.a_p", sbs_hcal_a_p, &b_sbs_hcal_a_p);
   fChain->SetBranchAddress("Ndata.sbs.hcal.a_time", &Ndata_sbs_hcal_a_time, &b_Ndata_sbs_hcal_a_time);
   fChain->SetBranchAddress("sbs.hcal.a_time", sbs_hcal_a_time, &b_sbs_hcal_a_time);
   fChain->SetBranchAddress("Ndata.sbs.hcal.adccol", &Ndata_sbs_hcal_adccol, &b_Ndata_sbs_hcal_adccol);
   fChain->SetBranchAddress("sbs.hcal.adccol", sbs_hcal_adccol, &b_sbs_hcal_adccol);
   fChain->SetBranchAddress("Ndata.sbs.hcal.adcelemID", &Ndata_sbs_hcal_adcelemID, &b_Ndata_sbs_hcal_adcelemID);
   fChain->SetBranchAddress("sbs.hcal.adcelemID", sbs_hcal_adcelemID, &b_sbs_hcal_adcelemID);
   fChain->SetBranchAddress("Ndata.sbs.hcal.adclayer", &Ndata_sbs_hcal_adclayer, &b_Ndata_sbs_hcal_adclayer);
   fChain->SetBranchAddress("sbs.hcal.adclayer", sbs_hcal_adclayer, &b_sbs_hcal_adclayer);
   fChain->SetBranchAddress("Ndata.sbs.hcal.adcrow", &Ndata_sbs_hcal_adcrow, &b_Ndata_sbs_hcal_adcrow);
   fChain->SetBranchAddress("sbs.hcal.adcrow", sbs_hcal_adcrow, &b_sbs_hcal_adcrow);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.col", &Ndata_sbs_hcal_clus_col, &b_Ndata_sbs_hcal_clus_col);
   fChain->SetBranchAddress("sbs.hcal.clus.col", sbs_hcal_clus_col, &b_sbs_hcal_clus_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.e", &Ndata_sbs_hcal_clus_e, &b_Ndata_sbs_hcal_clus_e);
   fChain->SetBranchAddress("sbs.hcal.clus.e", sbs_hcal_clus_e, &b_sbs_hcal_clus_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.e_c", &Ndata_sbs_hcal_clus_e_c, &b_Ndata_sbs_hcal_clus_e_c);
   fChain->SetBranchAddress("sbs.hcal.clus.e_c", sbs_hcal_clus_e_c, &b_sbs_hcal_clus_e_c);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.eblk", &Ndata_sbs_hcal_clus_eblk, &b_Ndata_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("sbs.hcal.clus.eblk", sbs_hcal_clus_eblk, &b_sbs_hcal_clus_eblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.eblk_c", &Ndata_sbs_hcal_clus_eblk_c, &b_Ndata_sbs_hcal_clus_eblk_c);
   fChain->SetBranchAddress("sbs.hcal.clus.eblk_c", sbs_hcal_clus_eblk_c, &b_sbs_hcal_clus_eblk_c);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id, &b_Ndata_sbs_hcal_clus_id);
   fChain->SetBranchAddress("sbs.hcal.clus.id", sbs_hcal_clus_id, &b_sbs_hcal_clus_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.nblk", &Ndata_sbs_hcal_clus_nblk, &b_Ndata_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("sbs.hcal.clus.nblk", sbs_hcal_clus_nblk, &b_sbs_hcal_clus_nblk);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.row", &Ndata_sbs_hcal_clus_row, &b_Ndata_sbs_hcal_clus_row);
   fChain->SetBranchAddress("sbs.hcal.clus.row", sbs_hcal_clus_row, &b_sbs_hcal_clus_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.x", &Ndata_sbs_hcal_clus_x, &b_Ndata_sbs_hcal_clus_x);
   fChain->SetBranchAddress("sbs.hcal.clus.x", sbs_hcal_clus_x, &b_sbs_hcal_clus_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus.y", &Ndata_sbs_hcal_clus_y, &b_Ndata_sbs_hcal_clus_y);
   fChain->SetBranchAddress("sbs.hcal.clus.y", sbs_hcal_clus_y, &b_sbs_hcal_clus_y);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.col", &Ndata_sbs_hcal_clus_blk_col, &b_Ndata_sbs_hcal_clus_blk_col);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.col", sbs_hcal_clus_blk_col, &b_sbs_hcal_clus_blk_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.e", &Ndata_sbs_hcal_clus_blk_e, &b_Ndata_sbs_hcal_clus_blk_e);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.e", sbs_hcal_clus_blk_e, &b_sbs_hcal_clus_blk_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.e_c", &Ndata_sbs_hcal_clus_blk_e_c, &b_Ndata_sbs_hcal_clus_blk_e_c);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.e_c", &sbs_hcal_clus_blk_e_c, &b_sbs_hcal_clus_blk_e_c);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.id", &Ndata_sbs_hcal_clus_blk_id, &b_Ndata_sbs_hcal_clus_blk_id);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.id", sbs_hcal_clus_blk_id, &b_sbs_hcal_clus_blk_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.row", &Ndata_sbs_hcal_clus_blk_row, &b_Ndata_sbs_hcal_clus_blk_row);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.row", sbs_hcal_clus_blk_row, &b_sbs_hcal_clus_blk_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.x", &Ndata_sbs_hcal_clus_blk_x, &b_Ndata_sbs_hcal_clus_blk_x);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.x", sbs_hcal_clus_blk_x, &b_sbs_hcal_clus_blk_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.clus_blk.y", &Ndata_sbs_hcal_clus_blk_y, &b_Ndata_sbs_hcal_clus_blk_y);
   fChain->SetBranchAddress("sbs.hcal.clus_blk.y", sbs_hcal_clus_blk_y, &b_sbs_hcal_clus_blk_y);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.atime", &Ndata_sbs_hcal_goodblock_atime, &b_Ndata_sbs_hcal_goodblock_atime);
   fChain->SetBranchAddress("sbs.hcal.goodblock.atime", sbs_hcal_goodblock_atime, &b_sbs_hcal_goodblock_atime);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.col", &Ndata_sbs_hcal_goodblock_col, &b_Ndata_sbs_hcal_goodblock_col);
   fChain->SetBranchAddress("sbs.hcal.goodblock.col", sbs_hcal_goodblock_col, &b_sbs_hcal_goodblock_col);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.e", &Ndata_sbs_hcal_goodblock_e, &b_Ndata_sbs_hcal_goodblock_e);
   fChain->SetBranchAddress("sbs.hcal.goodblock.e", sbs_hcal_goodblock_e, &b_sbs_hcal_goodblock_e);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.id", &Ndata_sbs_hcal_goodblock_id, &b_Ndata_sbs_hcal_goodblock_id);
   fChain->SetBranchAddress("sbs.hcal.goodblock.id", sbs_hcal_goodblock_id, &b_sbs_hcal_goodblock_id);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.row", &Ndata_sbs_hcal_goodblock_row, &b_Ndata_sbs_hcal_goodblock_row);
   fChain->SetBranchAddress("sbs.hcal.goodblock.row", sbs_hcal_goodblock_row, &b_sbs_hcal_goodblock_row);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.x", &Ndata_sbs_hcal_goodblock_x, &b_Ndata_sbs_hcal_goodblock_x);
   fChain->SetBranchAddress("sbs.hcal.goodblock.x", sbs_hcal_goodblock_x, &b_sbs_hcal_goodblock_x);
   fChain->SetBranchAddress("Ndata.sbs.hcal.goodblock.y", &Ndata_sbs_hcal_goodblock_y, &b_Ndata_sbs_hcal_goodblock_y);
   fChain->SetBranchAddress("sbs.hcal.goodblock.y", sbs_hcal_goodblock_y, &b_sbs_hcal_goodblock_y);
   fChain->SetBranchAddress("Ndata.sbs.hcal.ped", &Ndata_sbs_hcal_ped, &b_Ndata_sbs_hcal_ped);
   fChain->SetBranchAddress("sbs.hcal.ped", sbs_hcal_ped, &b_sbs_hcal_ped);
   fChain->SetBranchAddress("Ndata.sbs.hcal.tdc", &Ndata_sbs_hcal_tdc, &b_Ndata_sbs_hcal_tdc);
   fChain->SetBranchAddress("sbs.hcal.tdc", sbs_hcal_tdc, &b_sbs_hcal_tdc);
   fChain->SetBranchAddress("Ndata.sbs.hcal.tdc_mult", &Ndata_sbs_hcal_tdc_mult, &b_Ndata_sbs_hcal_tdc_mult);
   fChain->SetBranchAddress("sbs.hcal.tdc_mult", sbs_hcal_tdc_mult, &b_sbs_hcal_tdc_mult);
   fChain->SetBranchAddress("Ndata.sbs.hcal.tdccol", &Ndata_sbs_hcal_tdccol, &b_Ndata_sbs_hcal_tdccol);
   fChain->SetBranchAddress("sbs.hcal.tdccol", sbs_hcal_tdccol, &b_sbs_hcal_tdccol);
   fChain->SetBranchAddress("Ndata.sbs.hcal.tdcelemID", &Ndata_sbs_hcal_tdcelemID, &b_Ndata_sbs_hcal_tdcelemID);
   fChain->SetBranchAddress("sbs.hcal.tdcelemID", sbs_hcal_tdcelemID, &b_sbs_hcal_tdcelemID);
   fChain->SetBranchAddress("Ndata.sbs.hcal.tdclayer", &Ndata_sbs_hcal_tdclayer, &b_Ndata_sbs_hcal_tdclayer);
   fChain->SetBranchAddress("sbs.hcal.tdclayer", sbs_hcal_tdclayer, &b_sbs_hcal_tdclayer);
   fChain->SetBranchAddress("Ndata.sbs.hcal.tdcrow", &Ndata_sbs_hcal_tdcrow, &b_Ndata_sbs_hcal_tdcrow);
   fChain->SetBranchAddress("sbs.hcal.tdcrow", sbs_hcal_tdcrow, &b_sbs_hcal_tdcrow);
   fChain->SetBranchAddress("Ndata.sbs.tr.beta", &Ndata_sbs_tr_beta, &b_Ndata_sbs_tr_beta);
   fChain->SetBranchAddress("sbs.tr.beta", &sbs_tr_beta, &b_sbs_tr_beta);
   fChain->SetBranchAddress("Ndata.sbs.tr.chi2", &Ndata_sbs_tr_chi2, &b_Ndata_sbs_tr_chi2);
   fChain->SetBranchAddress("sbs.tr.chi2", &sbs_tr_chi2, &b_sbs_tr_chi2);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_ph", &Ndata_sbs_tr_d_ph, &b_Ndata_sbs_tr_d_ph);
   fChain->SetBranchAddress("sbs.tr.d_ph", &sbs_tr_d_ph, &b_sbs_tr_d_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_th", &Ndata_sbs_tr_d_th, &b_Ndata_sbs_tr_d_th);
   fChain->SetBranchAddress("sbs.tr.d_th", &sbs_tr_d_th, &b_sbs_tr_d_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_x", &Ndata_sbs_tr_d_x, &b_Ndata_sbs_tr_d_x);
   fChain->SetBranchAddress("sbs.tr.d_x", &sbs_tr_d_x, &b_sbs_tr_d_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.d_y", &Ndata_sbs_tr_d_y, &b_Ndata_sbs_tr_d_y);
   fChain->SetBranchAddress("sbs.tr.d_y", &sbs_tr_d_y, &b_sbs_tr_d_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.dbeta", &Ndata_sbs_tr_dbeta, &b_Ndata_sbs_tr_dbeta);
   fChain->SetBranchAddress("sbs.tr.dbeta", &sbs_tr_dbeta, &b_sbs_tr_dbeta);
   fChain->SetBranchAddress("Ndata.sbs.tr.dtime", &Ndata_sbs_tr_dtime, &b_Ndata_sbs_tr_dtime);
   fChain->SetBranchAddress("sbs.tr.dtime", &sbs_tr_dtime, &b_sbs_tr_dtime);
   fChain->SetBranchAddress("Ndata.sbs.tr.flag", &Ndata_sbs_tr_flag, &b_Ndata_sbs_tr_flag);
   fChain->SetBranchAddress("sbs.tr.flag", &sbs_tr_flag, &b_sbs_tr_flag);
   fChain->SetBranchAddress("Ndata.sbs.tr.ndof", &Ndata_sbs_tr_ndof, &b_Ndata_sbs_tr_ndof);
   fChain->SetBranchAddress("sbs.tr.ndof", &sbs_tr_ndof, &b_sbs_tr_ndof);
   fChain->SetBranchAddress("Ndata.sbs.tr.p", &Ndata_sbs_tr_p, &b_Ndata_sbs_tr_p);
   fChain->SetBranchAddress("sbs.tr.p", &sbs_tr_p, &b_sbs_tr_p);
   fChain->SetBranchAddress("Ndata.sbs.tr.pathl", &Ndata_sbs_tr_pathl, &b_Ndata_sbs_tr_pathl);
   fChain->SetBranchAddress("sbs.tr.pathl", &sbs_tr_pathl, &b_sbs_tr_pathl);
   fChain->SetBranchAddress("Ndata.sbs.tr.ph", &Ndata_sbs_tr_ph, &b_Ndata_sbs_tr_ph);
   fChain->SetBranchAddress("sbs.tr.ph", &sbs_tr_ph, &b_sbs_tr_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.px", &Ndata_sbs_tr_px, &b_Ndata_sbs_tr_px);
   fChain->SetBranchAddress("sbs.tr.px", &sbs_tr_px, &b_sbs_tr_px);
   fChain->SetBranchAddress("Ndata.sbs.tr.py", &Ndata_sbs_tr_py, &b_Ndata_sbs_tr_py);
   fChain->SetBranchAddress("sbs.tr.py", &sbs_tr_py, &b_sbs_tr_py);
   fChain->SetBranchAddress("Ndata.sbs.tr.pz", &Ndata_sbs_tr_pz, &b_Ndata_sbs_tr_pz);
   fChain->SetBranchAddress("sbs.tr.pz", &sbs_tr_pz, &b_sbs_tr_pz);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_ph", &Ndata_sbs_tr_r_ph, &b_Ndata_sbs_tr_r_ph);
   fChain->SetBranchAddress("sbs.tr.r_ph", &sbs_tr_r_ph, &b_sbs_tr_r_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_th", &Ndata_sbs_tr_r_th, &b_Ndata_sbs_tr_r_th);
   fChain->SetBranchAddress("sbs.tr.r_th", &sbs_tr_r_th, &b_sbs_tr_r_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_x", &Ndata_sbs_tr_r_x, &b_Ndata_sbs_tr_r_x);
   fChain->SetBranchAddress("sbs.tr.r_x", &sbs_tr_r_x, &b_sbs_tr_r_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.r_y", &Ndata_sbs_tr_r_y, &b_Ndata_sbs_tr_r_y);
   fChain->SetBranchAddress("sbs.tr.r_y", &sbs_tr_r_y, &b_sbs_tr_r_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_dp", &Ndata_sbs_tr_tg_dp, &b_Ndata_sbs_tr_tg_dp);
   fChain->SetBranchAddress("sbs.tr.tg_dp", &sbs_tr_tg_dp, &b_sbs_tr_tg_dp);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_ph", &Ndata_sbs_tr_tg_ph, &b_Ndata_sbs_tr_tg_ph);
   fChain->SetBranchAddress("sbs.tr.tg_ph", &sbs_tr_tg_ph, &b_sbs_tr_tg_ph);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_th", &Ndata_sbs_tr_tg_th, &b_Ndata_sbs_tr_tg_th);
   fChain->SetBranchAddress("sbs.tr.tg_th", &sbs_tr_tg_th, &b_sbs_tr_tg_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.tg_y", &Ndata_sbs_tr_tg_y, &b_Ndata_sbs_tr_tg_y);
   fChain->SetBranchAddress("sbs.tr.tg_y", &sbs_tr_tg_y, &b_sbs_tr_tg_y);
   fChain->SetBranchAddress("Ndata.sbs.tr.th", &Ndata_sbs_tr_th, &b_Ndata_sbs_tr_th);
   fChain->SetBranchAddress("sbs.tr.th", &sbs_tr_th, &b_sbs_tr_th);
   fChain->SetBranchAddress("Ndata.sbs.tr.time", &Ndata_sbs_tr_time, &b_Ndata_sbs_tr_time);
   fChain->SetBranchAddress("sbs.tr.time", &sbs_tr_time, &b_sbs_tr_time);
   fChain->SetBranchAddress("Ndata.sbs.tr.vx", &Ndata_sbs_tr_vx, &b_Ndata_sbs_tr_vx);
   fChain->SetBranchAddress("sbs.tr.vx", &sbs_tr_vx, &b_sbs_tr_vx);
   fChain->SetBranchAddress("Ndata.sbs.tr.vy", &Ndata_sbs_tr_vy, &b_Ndata_sbs_tr_vy);
   fChain->SetBranchAddress("sbs.tr.vy", &sbs_tr_vy, &b_sbs_tr_vy);
   fChain->SetBranchAddress("Ndata.sbs.tr.vz", &Ndata_sbs_tr_vz, &b_Ndata_sbs_tr_vz);
   fChain->SetBranchAddress("sbs.tr.vz", &sbs_tr_vz, &b_sbs_tr_vz);
   fChain->SetBranchAddress("Ndata.sbs.tr.x", &Ndata_sbs_tr_x, &b_Ndata_sbs_tr_x);
   fChain->SetBranchAddress("sbs.tr.x", &sbs_tr_x, &b_sbs_tr_x);
   fChain->SetBranchAddress("Ndata.sbs.tr.y", &Ndata_sbs_tr_y, &b_Ndata_sbs_tr_y);
   fChain->SetBranchAddress("sbs.tr.y", &sbs_tr_y, &b_sbs_tr_y);
   fChain->SetBranchAddress("MC.hit.n", &MC_hit_n, &b_MC_hit_n);
   fChain->SetBranchAddress("MC.mc_px", &MC_mc_px, &b_MC_mc_px);
   fChain->SetBranchAddress("MC.mc_py", &MC_mc_py, &b_MC_mc_py);
   fChain->SetBranchAddress("MC.mc_pz", &MC_mc_pz, &b_MC_mc_pz);
   fChain->SetBranchAddress("MC.mc_vx", &MC_mc_vx, &b_MC_mc_vx);
   fChain->SetBranchAddress("MC.mc_vy", &MC_mc_vy, &b_MC_mc_vy);
   fChain->SetBranchAddress("MC.mc_vz", &MC_mc_vz, &b_MC_mc_vz);
   fChain->SetBranchAddress("MC.nbbgemhits", &MC_nbbgemhits, &b_MC_nbbgemhits);
   fChain->SetBranchAddress("MC.nbbtracks", &MC_nbbtracks, &b_MC_nbbtracks);
   fChain->SetBranchAddress("MC.pt.n", &MC_pt_n, &b_MC_pt_n);
   fChain->SetBranchAddress("MC.tr.n", &MC_tr_n, &b_MC_tr_n);
   fChain->SetBranchAddress("MC.weight", &MC_weight, &b_MC_weight);
   fChain->SetBranchAddress("bb.gem.hit.ngoodhits", &bb_gem_hit_ngoodhits, &b_bb_gem_hit_ngoodhits);
   fChain->SetBranchAddress("bb.gem.nlayershit", &bb_gem_nlayershit, &b_bb_gem_nlayershit);
   fChain->SetBranchAddress("bb.gem.nlayershitu", &bb_gem_nlayershitu, &b_bb_gem_nlayershitu);
   fChain->SetBranchAddress("bb.gem.nlayershituv", &bb_gem_nlayershituv, &b_bb_gem_nlayershituv);
   fChain->SetBranchAddress("bb.gem.nlayershitv", &bb_gem_nlayershitv, &b_bb_gem_nlayershitv);
   fChain->SetBranchAddress("bb.gem.track.besttrack", &bb_gem_track_besttrack, &b_bb_gem_track_besttrack);
   fChain->SetBranchAddress("bb.gem.track.ntrack", &bb_gem_track_ntrack, &b_bb_gem_track_ntrack);
   fChain->SetBranchAddress("bb.grinch.nclus", &bb_grinch_nclus, &b_bb_grinch_nclus);
   fChain->SetBranchAddress("bb.grinch.nhits", &bb_grinch_nhits, &b_bb_grinch_nhits);
   fChain->SetBranchAddress("bb.hodo.bar.ngoodbars", &bb_hodo_bar_ngoodbars, &b_bb_hodo_bar_ngoodbars);
   fChain->SetBranchAddress("bb.hodo.nclus", &bb_hodo_nclus, &b_bb_hodo_nclus);
   fChain->SetBranchAddress("bb.hodo.ngoodADChits", &bb_hodo_ngoodADChits, &b_bb_hodo_ngoodADChits);
   fChain->SetBranchAddress("bb.hodo.ngoodTDChits", &bb_hodo_ngoodTDChits, &b_bb_hodo_ngoodTDChits);
   fChain->SetBranchAddress("bb.hodo.nhits", &bb_hodo_nhits, &b_bb_hodo_nhits);
   fChain->SetBranchAddress("bb.hodo.nrefhits", &bb_hodo_nrefhits, &b_bb_hodo_nrefhits);
   fChain->SetBranchAddress("bb.ps.colblk", &bb_ps_colblk, &b_bb_ps_colblk);
   fChain->SetBranchAddress("bb.ps.e", &bb_ps_e, &b_bb_ps_e);
   fChain->SetBranchAddress("bb.ps.e_c", &bb_ps_e_c, &b_bb_ps_e_c);
   fChain->SetBranchAddress("bb.ps.e_m_res", &bb_ps_e_m_res, &b_bb_ps_e_m_res);
   fChain->SetBranchAddress("bb.ps.eblk", &bb_ps_eblk, &b_bb_ps_eblk);
   fChain->SetBranchAddress("bb.ps.eblk_c", &bb_ps_eblk_c, &b_bb_ps_eblk_c);
   fChain->SetBranchAddress("bb.ps.idblk", &bb_ps_idblk, &b_bb_ps_idblk);
   fChain->SetBranchAddress("bb.ps.nblk", &bb_ps_nblk, &b_bb_ps_nblk);
   fChain->SetBranchAddress("bb.ps.nclus", &bb_ps_nclus, &b_bb_ps_nclus);
   fChain->SetBranchAddress("bb.ps.ngoodADChits", &bb_ps_ngoodADChits, &b_bb_ps_ngoodADChits);
   fChain->SetBranchAddress("bb.ps.ngoodTDChits", &bb_ps_ngoodTDChits, &b_bb_ps_ngoodTDChits);
   fChain->SetBranchAddress("bb.ps.nhits", &bb_ps_nhits, &b_bb_ps_nhits);
   fChain->SetBranchAddress("bb.ps.nrefhits", &bb_ps_nrefhits, &b_bb_ps_nrefhits);
   fChain->SetBranchAddress("bb.ps.rowblk", &bb_ps_rowblk, &b_bb_ps_rowblk);
   fChain->SetBranchAddress("bb.ps.x", &bb_ps_x, &b_bb_ps_x);
   fChain->SetBranchAddress("bb.ps.x_m_res", &bb_ps_x_m_res, &b_bb_ps_x_m_res);
   fChain->SetBranchAddress("bb.ps.y", &bb_ps_y, &b_bb_ps_y);
   fChain->SetBranchAddress("bb.ps.y_m_res", &bb_ps_y_m_res, &b_bb_ps_y_m_res);
   fChain->SetBranchAddress("bb.sh.colblk", &bb_sh_colblk, &b_bb_sh_colblk);
   fChain->SetBranchAddress("bb.sh.e", &bb_sh_e, &b_bb_sh_e);
   fChain->SetBranchAddress("bb.sh.e_c", &bb_sh_e_c, &b_bb_sh_e_c);
   fChain->SetBranchAddress("bb.sh.e_m_res", &bb_sh_e_m_res, &b_bb_sh_e_m_res);
   fChain->SetBranchAddress("bb.sh.eblk", &bb_sh_eblk, &b_bb_sh_eblk);
   fChain->SetBranchAddress("bb.sh.eblk_c", &bb_sh_eblk_c, &b_bb_sh_eblk_c);
   fChain->SetBranchAddress("bb.sh.idblk", &bb_sh_idblk, &b_bb_sh_idblk);
   fChain->SetBranchAddress("bb.sh.nblk", &bb_sh_nblk, &b_bb_sh_nblk);
   fChain->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus, &b_bb_sh_nclus);
   fChain->SetBranchAddress("bb.sh.ngoodADChits", &bb_sh_ngoodADChits, &b_bb_sh_ngoodADChits);
   fChain->SetBranchAddress("bb.sh.ngoodTDChits", &bb_sh_ngoodTDChits, &b_bb_sh_ngoodTDChits);
   fChain->SetBranchAddress("bb.sh.nhits", &bb_sh_nhits, &b_bb_sh_nhits);
   fChain->SetBranchAddress("bb.sh.nrefhits", &bb_sh_nrefhits, &b_bb_sh_nrefhits);
   fChain->SetBranchAddress("bb.sh.rowblk", &bb_sh_rowblk, &b_bb_sh_rowblk);
   fChain->SetBranchAddress("bb.sh.x", &bb_sh_x, &b_bb_sh_x);
   fChain->SetBranchAddress("bb.sh.x_m_res", &bb_sh_x_m_res, &b_bb_sh_x_m_res);
   fChain->SetBranchAddress("bb.sh.y", &bb_sh_y, &b_bb_sh_y);
   fChain->SetBranchAddress("bb.sh.y_m_res", &bb_sh_y_m_res, &b_bb_sh_y_m_res);
   fChain->SetBranchAddress("bb.tr.n", &bb_tr_n, &b_bb_tr_n);
   fChain->SetBranchAddress("sbs.hcal.colblk", &sbs_hcal_colblk, &b_sbs_hcal_colblk);
   fChain->SetBranchAddress("sbs.hcal.e", &sbs_hcal_e, &b_sbs_hcal_e);
   fChain->SetBranchAddress("sbs.hcal.e_c", &sbs_hcal_e_c, &b_sbs_hcal_e_c);
   fChain->SetBranchAddress("sbs.hcal.eblk", &sbs_hcal_eblk, &b_sbs_hcal_eblk);
   fChain->SetBranchAddress("sbs.hcal.eblk_c", &sbs_hcal_eblk_c, &b_sbs_hcal_eblk_c);
   fChain->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_idblk, &b_sbs_hcal_idblk);
   fChain->SetBranchAddress("sbs.hcal.ledbit", &sbs_hcal_ledbit, &b_sbs_hcal_ledbit);
   fChain->SetBranchAddress("sbs.hcal.ledcount", &sbs_hcal_ledcount, &b_sbs_hcal_ledcount);
   fChain->SetBranchAddress("sbs.hcal.nblk", &sbs_hcal_nblk, &b_sbs_hcal_nblk);
   fChain->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus, &b_sbs_hcal_nclus);
   fChain->SetBranchAddress("sbs.hcal.ngoodADChits", &sbs_hcal_ngoodADChits, &b_sbs_hcal_ngoodADChits);
   fChain->SetBranchAddress("sbs.hcal.ngoodTDChits", &sbs_hcal_ngoodTDChits, &b_sbs_hcal_ngoodTDChits);
   fChain->SetBranchAddress("sbs.hcal.nhits", &sbs_hcal_nhits, &b_sbs_hcal_nhits);
   fChain->SetBranchAddress("sbs.hcal.nrefhits", &sbs_hcal_nrefhits, &b_sbs_hcal_nrefhits);
   fChain->SetBranchAddress("sbs.hcal.rowblk", &sbs_hcal_rowblk, &b_sbs_hcal_rowblk);
   fChain->SetBranchAddress("sbs.hcal.x", &sbs_hcal_x, &b_sbs_hcal_x);
   fChain->SetBranchAddress("sbs.hcal.y", &sbs_hcal_y, &b_sbs_hcal_y);
   fChain->SetBranchAddress("sbs.status", &sbs_status, &b_sbs_status);
   fChain->SetBranchAddress("sbs.tr.n", &sbs_tr_n, &b_sbs_tr_n);
   fChain->SetBranchAddress("fEvtHdr.fEvtTime", &fEvtHdr_fEvtTime, &b_Event_Branch_fEvtHdr_fEvtTime);
   fChain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtHdr_fEvtNum, &b_Event_Branch_fEvtHdr_fEvtNum);
   fChain->SetBranchAddress("fEvtHdr.fEvtType", &fEvtHdr_fEvtType, &b_Event_Branch_fEvtHdr_fEvtType);
   fChain->SetBranchAddress("fEvtHdr.fEvtLen", &fEvtHdr_fEvtLen, &b_Event_Branch_fEvtHdr_fEvtLen);
   fChain->SetBranchAddress("fEvtHdr.fHelicity", &fEvtHdr_fHelicity, &b_Event_Branch_fEvtHdr_fHelicity);
   fChain->SetBranchAddress("fEvtHdr.fTargetPol", &fEvtHdr_fTargetPol, &b_Event_Branch_fEvtHdr_fTargetPol);
   fChain->SetBranchAddress("fEvtHdr.fRun", &fEvtHdr_fRun, &b_Event_Branch_fEvtHdr_fRun);
   Notify();
}

Bool_t gmn_rec_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gmn_rec_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gmn_rec_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gmn_rec_tree_cxx
