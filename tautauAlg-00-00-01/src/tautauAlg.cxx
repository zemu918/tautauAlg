//-----------------------------
// e+ e- -> tau+ tau-,
//          |    |-> e- nu anti-nu
//          mu+ nu anti_nu
//----------------------------

#include "tautauAlg/util/MyConst.h"
#include "tautauAlg/util/MyUtil.h"
#include "tautauAlg/util/MyIsGoodshower.h"
#include "tautauAlg/util/MyIsGoodtrack.h"
#include "tautauAlg/util/MyInitIP.h"
#include "tautauAlg/util/MyInfoShower.h"
#include "tautauAlg/util/MyInfoCharge.h"
#include "tautauAlg/util/MyGamConv.h"
#include "tautauAlg/util/Mygamee.h"
#include "tautauAlg/tautauAlg.h"

using namespace Event;

typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
int Ncut[10];

//***************************Declare******************************
tautauAlg::tautauAlg(const std::string &name, ISvcLocator *pSvcLocator)
    : Algorithm(name, pSvcLocator) {
  m_nEvtDisp = 100;
  m_nEvt = 0;
  declareProperty("writeGenOnly", m_writeGenOnly = false);
  declareProperty("Ecms", m_ecms = 4.26);
  declareProperty("FSRCor", m_appFSRCorrection = true);
  declareProperty("m_testMC", m_testMC = 1);
}

//***************************Fill Tree******************************
StatusCode tautauAlg::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in initialize()" << endmsg;

  StatusCode status;

  NTuplePtr nt1(ntupleSvc(), "FILE1/mctruth");
  if (nt1)
    m_tuple1 = nt1;
  else {
    m_tuple1 =
        ntupleSvc()->book("FILE1/mctruth", CLID_ColumnWiseTuple, "mctruth");
    if (m_tuple1) {
      m_tuple1->addItem("run", m_run);
      m_tuple1->addItem("event", m_event);
      m_tuple1->addItem("indexmc", m_idxmc, 0, 100);
      m_tuple1->addIndexedItem("pdgid", m_idxmc, m_pdgid);
      m_tuple1->addIndexedItem("motheridx", m_idxmc, m_motheridx);
      m_tuple1->addItem("info_tru", 7, 8, m_info_tru);
    } else {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1)
          << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt3(ntupleSvc(), "FILE1/gen");
  if (nt3)
    m_tuple3 = nt3;
  else {
    m_tuple3 =
        ntupleSvc()->book("FILE1/gen", CLID_ColumnWiseTuple, "gen");
    if (m_tuple3) {
      m_tuple3->addItem("within_kine_bounds", m_within_kine_bounds); 
      m_tuple3->addItem("pion_px", m_pi_px);
      m_tuple3->addItem("pion_py", m_pi_py);
      m_tuple3->addItem("pion_pz", m_pi_pz);
      m_tuple3->addItem("pion_e", m_pi_e);

      m_tuple3->addItem("gam_px", m_gam_px);
      m_tuple3->addItem("gam_py", m_gam_py);
      m_tuple3->addItem("gam_pz", m_gam_pz);
      m_tuple3->addItem("gam_e", m_gam_e);

      m_tuple3->addItem("tau_px", m_tau_px);
      m_tuple3->addItem("tau_py", m_tau_py);
      m_tuple3->addItem("tau_pz", m_tau_pz);
      m_tuple3->addItem("tau_e", m_tau_e);

      m_tuple3->addItem("nu_px", m_nu_px);
      m_tuple3->addItem("nu_py", m_nu_py);
      m_tuple3->addItem("nu_pz", m_nu_pz);
      m_tuple3->addItem("nu_e", m_nu_e);

      m_tuple3->addItem("x", m_x);
      m_tuple3->addItem("y", m_y);
      m_tuple3->addItem("t", m_t);
      m_tuple3->addItem("s", m_s);
      m_tuple3->addItem("E", m_E);
      m_tuple3->addItem("z", m_z);

      m_tuple3->addItem("f_IB", m_f_IB);
      m_tuple3->addItem("f_VV", m_f_VV);
      m_tuple3->addItem("f_AV", m_f_AV);
      m_tuple3->addItem("f_AA", m_f_AA);
      m_tuple3->addItem("f_IB_V", m_f_IB_V);
      m_tuple3->addItem("f_IB_A", m_f_IB_A);
      m_tuple3->addItem("Gam_IB", m_Gam_IB);
      m_tuple3->addItem("Gam_SD", m_Gam_SD);
      m_tuple3->addItem("Gam_int", m_Gam_int);
    } else {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple3)
          << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt2(ntupleSvc(), "FILE1/cccc");
  if (nt2)
    m_tuple2 = nt2;
  else {
    m_tuple2 = ntupleSvc()->book("FILE1/cccc", CLID_ColumnWiseTuple,
                                 " N-Tuple example");
    if (m_tuple2) {
      m_tuple2->addItem("rec_run", m_rec_run);
      m_tuple2->addItem("rec_event", m_rec_event);
      m_tuple2->addItem("rec_nchrg", m_rec_nchrg);
      m_tuple2->addItem("rec_nneu", m_rec_nneu);

      m_tuple2->addItem("indexmc", m_rec_idxmc, 0, 100);
      m_tuple2->addIndexedItem("pdgid", m_rec_idxmc, m_rec_pdgid);
      m_tuple2->addIndexedItem("motheridx", m_rec_idxmc, m_rec_motheridx);

      // 2020.03.14-15,
      //*****************
      // MDC information
      m_tuple2->addItem("charge_1", m_charge_1);
      m_tuple2->addItem("charge_2", m_charge_2);
      m_tuple2->addItem("pmdc_1", m_pmdc_1);
      m_tuple2->addItem("pmdc_2", m_pmdc_2);
      m_tuple2->addItem("ptmdc_1", m_ptmdc_1);
      m_tuple2->addItem("ptmdc_2", m_ptmdc_2);
      m_tuple2->addItem("pkal_1", 5, m_pkal_1);
      m_tuple2->addItem("pkal_2", 5, m_pkal_2);
      m_tuple2->addItem("ptkal_1", 5, m_ptkal_1);
      m_tuple2->addItem("ptkal_2", 5, m_ptkal_2);
      m_tuple2->addItem("p3_1", 3, m_p3_1);
      m_tuple2->addItem("p3_2", 3, m_p3_2);
      m_tuple2->addItem("theta_1", m_theta_1);
      m_tuple2->addItem("theta_2", m_theta_2);

      // PID information
      m_tuple2->addItem("probe_1", m_probe_1);
      m_tuple2->addItem("probmu_1", m_probmu_1);
      m_tuple2->addItem("probpi_1", m_probpi_1);
      m_tuple2->addItem("probk_1", m_probk_1);
      m_tuple2->addItem("probp_1", m_probp_1);
      m_tuple2->addItem("probe_2", m_probe_2);
      m_tuple2->addItem("probmu_2", m_probmu_2);
      m_tuple2->addItem("probpi_2", m_probpi_2);
      m_tuple2->addItem("probk_2", m_probk_2);
      m_tuple2->addItem("probp_2", m_probp_2);

      // TOF information
      m_tuple2->addItem("fgtof_1", m_fgtof_1);
      m_tuple2->addItem("counter_1", m_counter_1);
      m_tuple2->addItem("isbarrel_1", m_isbarrel_1);
      m_tuple2->addItem("layertof_1", m_layertof_1);
      m_tuple2->addItem("iscluster_1", m_iscluster_1);
      m_tuple2->addItem("tof_1", m_tof_1);
      m_tuple2->addItem("texe_1", m_texe_1);
      m_tuple2->addItem("texmu_1", m_texmu_1);
      m_tuple2->addItem("texpi_1", m_texpi_1);
      m_tuple2->addItem("texk_1", m_texk_1);
      m_tuple2->addItem("texp_1", m_texp_1);
      m_tuple2->addItem("dte_1", m_dte_1);
      m_tuple2->addItem("dtmu_1", m_dtmu_1);
      m_tuple2->addItem("dtpi_1", m_dtpi_1);
      m_tuple2->addItem("dtk_1", m_dtk_1);
      m_tuple2->addItem("dtp_1", m_dtp_1);
      m_tuple2->addItem("fgtof_2", m_fgtof_2);
      m_tuple2->addItem("counter_2", m_counter_2);
      m_tuple2->addItem("isbarrel_2", m_isbarrel_2);
      m_tuple2->addItem("layertof_2", m_layertof_2);
      m_tuple2->addItem("iscluster_2", m_iscluster_2);
      m_tuple2->addItem("tof_2", m_tof_2);
      m_tuple2->addItem("texe_2", m_texe_2);
      m_tuple2->addItem("texmu_2", m_texmu_2);
      m_tuple2->addItem("texpi_2", m_texpi_2);
      m_tuple2->addItem("texk_2", m_texk_2);
      m_tuple2->addItem("texp_2", m_texp_2);
      m_tuple2->addItem("dte_2", m_dte_2);
      m_tuple2->addItem("dtmu_2", m_dtmu_2);
      m_tuple2->addItem("dtpi_2", m_dtpi_2);

      m_tuple2->addItem("dtk_2", m_dtk_2);
      m_tuple2->addItem("dtp_2", m_dtp_2);

      // MUC information
      m_tuple2->addItem("maxhitsinlay_1", m_maxhitsinlay_1);
      m_tuple2->addItem("numhits_1", m_numhits_1);
      m_tuple2->addItem("numlayers_1", m_numlayers_1);
      m_tuple2->addItem("depth_1", m_depth_1);
      m_tuple2->addItem("mucchi2_1", m_mucchi2_1);
      m_tuple2->addItem("maxhitsinlay_2", m_maxhitsinlay_2);
      m_tuple2->addItem("numhits_2", m_numhits_2);
      m_tuple2->addItem("numlayers_2", m_numlayers_2);
      m_tuple2->addItem("depth_2", m_depth_2);
      m_tuple2->addItem("mucchi2_2", m_mucchi2_2);

      // EMC information
      m_tuple2->addItem("evp_1", m_evp_1);
      m_tuple2->addItem("ene_1", m_ene_1);
      m_tuple2->addItem("evp_2", m_evp_2);
      m_tuple2->addItem("ene_2", m_ene_2);
      m_tuple2->addItem("nTrk_EMC", m_nTrk_EMC);

      // dedx information
      m_tuple2->addItem("chi_e_1", m_chi_e_1);
      m_tuple2->addItem("chi_mu_1", m_chi_mu_1);
      m_tuple2->addItem("chi_pi_1", m_chi_pi_1);
      m_tuple2->addItem("chi_k_1", m_chi_k_1);
      m_tuple2->addItem("chi_p_1", m_chi_p_1);
      m_tuple2->addItem("chi_e_2", m_chi_e_2);
      m_tuple2->addItem("chi_mu_2", m_chi_mu_2);
      m_tuple2->addItem("chi_pi_2", m_chi_pi_2);
      m_tuple2->addItem("chi_k_2", m_chi_k_2);
      m_tuple2->addItem("chi_p_2", m_chi_p_2);

      m_tuple2->addItem("mTrks", 5, m_mTrks);
      m_tuple2->addItem("m2Trks", 5, m_m2Trks);
      m_tuple2->addItem("pkal_tot_px", 5, m_pkal_tot_px);
      m_tuple2->addItem("pkal_tot_py", 5, m_pkal_tot_py);
      m_tuple2->addItem("pkal_tot_pz", 5, m_pkal_tot_pz);
      m_tuple2->addItem("pkal_tot_e", 5, m_pkal_tot_e);
      m_tuple2->addItem("ang_mdc_acol", m_ang_mdc_acol);
      m_tuple2->addItem("ang_mdc_acop", m_ang_mdc_acop);
      m_tuple2->addItem("the_add", m_the_add);
      m_tuple2->addItem("phi_diff", m_phi_diff);
      m_tuple2->addItem("ptem", m_PTEM);
      m_tuple2->addItem("ee_acop", m_ee_acop);
      m_tuple2->addItem("emu_acop", m_emu_acop);
      m_tuple2->addItem("mue_acop", m_mue_acop);
      m_tuple2->addItem("epi_acop", m_epi_acop);
      m_tuple2->addItem("pie_acop", m_pie_acop);
      m_tuple2->addItem("ek_acop", m_ek_acop);
      m_tuple2->addItem("ke_acop", m_ke_acop);
      m_tuple2->addItem("mumu_acop", m_mumu_acop);
      m_tuple2->addItem("mupi_acop", m_mupi_acop);
      m_tuple2->addItem("pimu_acop", m_pimu_acop);
      m_tuple2->addItem("pipi_acop", m_pipi_acop);
      m_tuple2->addItem("muk_acop", m_muk_acop);
      m_tuple2->addItem("kmu_acop", m_kmu_acop);
      m_tuple2->addItem("pik_acop", m_pik_acop);
      m_tuple2->addItem("kpi_acop", m_kpi_acop);
      m_tuple2->addItem("kk_acop", m_kk_acop);
      m_tuple2->addItem("erho_11_acop", m_erho_11_acop);
      m_tuple2->addItem("erho_12_acop", m_erho_12_acop);
      m_tuple2->addItem("murho_11_acop", m_murho_11_acop);
      m_tuple2->addItem("murho_12_acop", m_murho_12_acop);
      m_tuple2->addItem("rhorho_11_acop", m_rhorho_11_acop);
      m_tuple2->addItem("rhorho_12_acop", m_rhorho_12_acop);
      m_tuple2->addItem("ee_PTEM", m_ee_PTEM);
      m_tuple2->addItem("emu_PTEM", m_emu_PTEM);
      m_tuple2->addItem("mue_PTEM", m_mue_PTEM);
      m_tuple2->addItem("epi_PTEM", m_epi_PTEM);
      m_tuple2->addItem("pie_PTEM", m_pie_PTEM);
      m_tuple2->addItem("ek_PTEM", m_ek_PTEM);
      m_tuple2->addItem("ke_PTEM", m_ke_PTEM);
      m_tuple2->addItem("mumu_PTEM", m_mumu_PTEM);
      m_tuple2->addItem("mupi_PTEM", m_mupi_PTEM);
      m_tuple2->addItem("pimu_PTEM", m_pimu_PTEM);
      m_tuple2->addItem("pipi_PTEM", m_pipi_PTEM);
      m_tuple2->addItem("muk_PTEM", m_muk_PTEM);

      m_tuple2->addItem("kmu_PTEM", m_kmu_PTEM);
      m_tuple2->addItem("pik_PTEM", m_pik_PTEM);
      m_tuple2->addItem("kpi_PTEM", m_kpi_PTEM);
      m_tuple2->addItem("kk_PTEM", m_kk_PTEM);
      m_tuple2->addItem("dlt_mpi0", m_dlt_mpi0);
      m_tuple2->addItem("eg_11", m_eg_11);
      m_tuple2->addItem("eg_12", m_eg_12);
      m_tuple2->addItem("eg_21", m_eg_21);
      m_tuple2->addItem("eg_22", m_eg_22);
      m_tuple2->addItem("mpi0_1", m_mpi0_1);
      m_tuple2->addItem("mpi0_2", m_mpi0_2);
      m_tuple2->addItem("mrho_11", m_mrho_11);
      m_tuple2->addItem("mrho_12", m_mrho_12);
      m_tuple2->addItem("mrho_21", m_mrho_21);
      m_tuple2->addItem("mrho_22", m_mrho_22);
      m_tuple2->addItem("prho_11", m_prho_11);
      m_tuple2->addItem("prho_12", m_prho_12);
      m_tuple2->addItem("prho_21", m_prho_21);
      m_tuple2->addItem("prho_22", m_prho_22);
      m_tuple2->addItem("theta_pi_11", m_theta_pi_11);
      m_tuple2->addItem("theta_pi_12", m_theta_pi_12);
      m_tuple2->addItem("theta_pi_21", m_theta_pi_21);
      m_tuple2->addItem("theta_pi_22", m_theta_pi_22);
      m_tuple2->addItem("ppi_cms_11", m_ppi_cms_11);
      m_tuple2->addItem("ppi_cms_12", m_ppi_cms_12);
      m_tuple2->addItem("ppi_cms_21", m_ppi_cms_21);
      m_tuple2->addItem("ppi_cms_22", m_ppi_cms_22);
      m_tuple2->addItem("theta_m_11", m_theta_m_11);
      m_tuple2->addItem("theta_m_12", m_theta_m_12);
      m_tuple2->addItem("theta_m_21", m_theta_m_21);
      m_tuple2->addItem("theta_m_22", m_theta_m_22);
      //*****************

      // 2020.04.02 after group meeting
      m_tuple2->addItem("ee_angle", m_ee_angle);
      m_tuple2->addItem("emu_angle", m_emu_angle);
      m_tuple2->addItem("mue_angle", m_mue_angle);
      m_tuple2->addItem("epi_angle", m_epi_angle);
      m_tuple2->addItem("pie_angle", m_pie_angle);
      m_tuple2->addItem("ek_angle", m_ek_angle);
      m_tuple2->addItem("ke_angle", m_ke_angle);
      m_tuple2->addItem("mumu_angle", m_mumu_angle);
      m_tuple2->addItem("mupi_angle", m_mupi_angle);
      m_tuple2->addItem("pimu_angle", m_pimu_angle);
      m_tuple2->addItem("pipi_angle", m_pipi_angle);
      m_tuple2->addItem("muk_angle", m_muk_angle);
      m_tuple2->addItem("kmu_angle", m_kmu_angle);
      m_tuple2->addItem("pik_angle", m_pik_angle);
      m_tuple2->addItem("kpi_angle", m_kpi_angle);
      m_tuple2->addItem("kk_angle", m_kk_angle);
      m_tuple2->addItem("erho_11_angle", m_erho_11_angle);
      m_tuple2->addItem("erho_12_angle", m_erho_12_angle);
      m_tuple2->addItem("murho_11_angle", m_murho_11_angle);
      m_tuple2->addItem("murho_12_angle", m_murho_12_angle);
      m_tuple2->addItem("rhorho_11_angle", m_rhorho_11_angle);
      m_tuple2->addItem("rhorho_12_angle", m_rhorho_12_angle);

      m_tuple2->addItem("miss_m2_ee", m_miss_m2_ee);
      m_tuple1->addItem("event", m_event);
      m_tuple2->addItem("miss_m2_emu", m_miss_m2_emu);
      m_tuple2->addItem("miss_m2_mue", m_miss_m2_mue);
      m_tuple2->addItem("miss_m2_epi", m_miss_m2_epi);
      m_tuple2->addItem("miss_m2_pie", m_miss_m2_pie);
      m_tuple2->addItem("miss_m2_ek", m_miss_m2_ek);
      m_tuple2->addItem("miss_m2_ke", m_miss_m2_ke);
      m_tuple2->addItem("miss_m2_mumu", m_miss_m2_mumu);
      m_tuple2->addItem("miss_m2_mupi", m_miss_m2_mupi);
      m_tuple2->addItem("miss_m2_pimu", m_miss_m2_pimu);
      m_tuple2->addItem("miss_m2_pipi", m_miss_m2_pipi);
      m_tuple2->addItem("miss_m2_muk", m_miss_m2_muk);
      m_tuple2->addItem("miss_m2_kmu", m_miss_m2_kmu);
      m_tuple2->addItem("miss_m2_pik", m_miss_m2_pik);
      m_tuple2->addItem("miss_m2_kpi", m_miss_m2_kpi);
      //    ma_nu_e = 0.0;
      m_tuple2->addItem("miss_m2_kk", m_miss_m2_kk);
      m_tuple2->addItem("miss_m2_erho_11", m_miss_m2_erho_11);
      //    ma_nu_e = 0.0;
      m_tuple2->addItem("miss_m2_erho_12", m_miss_m2_erho_12);
      m_tuple2->addItem("miss_m2_murho_11", m_miss_m2_murho_11);
      m_tuple2->addItem("miss_m2_murho_12", m_miss_m2_murho_12);
      m_tuple2->addItem("miss_m2_rhorho_11", m_miss_m2_rhorho_11);
      m_tuple2->addItem("miss_m2_rhorho_12", m_miss_m2_rhorho_12);

      m_tuple2->addItem("nsig_shower", m_nsig_shower, 0, 100);
      m_tuple2->addIndexedItem("shower_info", m_nsig_shower, 19, m_shower_info);
      m_tuple2->addItem("tot_egam", m_esum);
      m_tuple2->addItem("max_egam", m_egam_max);
      m_tuple2->addItem("max_gam_id", m_id_gam_max);
      m_tuple2->addItem("gamcon", 26, m_gamcon);
      m_tuple2->addItem("info_elp", 26, m_elp_info);
      m_tuple2->addItem("info_mum", 26, m_mum_info);
      m_tuple2->addItem("info_eta_gamee", 16, m_eta_gamee);
      m_tuple2->addItem("info_pi0_gamee", 16, m_pi0_gamee);
      m_tuple2->addItem("prob_elp", 6, m_elp_prob);
      m_tuple2->addItem("prob_mum", 6, m_mum_prob);

      m_tuple2->addItem("nelpFSR", m_nleppFSR, 0, 100);
      m_tuple2->addItem("elpFSRArray", m_nleppFSR, m_eleppFSRArray);
      m_tuple2->addItem("dtheelpFSRArray", m_nleppFSR, m_dtheleppFSRArray);

      m_tuple2->addItem("info_gam_fromelp", 8, m_gam_fromelp);
      m_tuple2->addItem("info_elp_corfsr", 8, m_elp_corfsr);
      m_tuple2->addItem("info_elp_mum", 3, m_info_elp_mum);
      //    ma_nu_e = 0.0;
      m_tuple2->addItem("info_p4_emu", 3, 8, m_info_emu);

      // 2020.03.27 for topology
      m_tuple2->addItem("nTrackMC", m_nTrack, 0, 100);
      m_tuple2->addIndexedItem("trackIDMC", m_nTrack, ma_trackID);
      m_tuple2->addIndexedItem("trackIndexMC", m_nTrack, ma_trackIndex);
      m_tuple2->addIndexedItem("motherIDMC", m_nTrack, ma_motherID);
      m_tuple2->addIndexedItem("motherIndexMC", m_nTrack, ma_motherIndex);
      m_tuple2->addIndexedItem("fromGeneratorMC", m_nTrack, ma_fromGenerator);
      m_tuple2->addIndexedItem("primaryParticleMC", m_nTrack,
                               ma_primaryParticle);
      m_tuple2->addIndexedItem("pxMC", m_nTrack, ma_px);
      m_tuple2->addIndexedItem("pyMC", m_nTrack, ma_py);
      m_tuple2->addIndexedItem("pzMC", m_nTrack, ma_pz);
      m_tuple2->addIndexedItem("eMC", m_nTrack, ma_e);

    } else {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple2)
          << endmsg;
      return StatusCode::FAILURE;
    }
  }
  log << MSG::INFO << "successfully return from initialize()" << endmsg;
  return StatusCode::SUCCESS;
}

//***************************Execution Event
//Filtering******************************
StatusCode tautauAlg::execute() {
  MsgStream log(msgSvc(), name());
  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),
                                               "/Event/EventHeader");
  if (!eventHeader) {
    log << MSG::FATAL << "Could not find Event Header" << endreq;

    return (StatusCode::FAILURE);
  }
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(),
                                        "/Event/EvtRec/EvtRecEvent");
  if (!evtRecEvent) {
    log << MSG::FATAL << "Could not find EvtRecEvent" << endreq;
    return StatusCode::FAILURE;
  }
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),
                                            "/Event/EvtRec/EvtRecTrackCol");
  if (!evtRecTrkCol) {
    log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
    return StatusCode::FAILURE;
  }

  EvtRecTrackIterator track_begin = evtRecTrkCol->begin();
  EvtRecTrackIterator charged_begin = track_begin;
  EvtRecTrackIterator charged_end = track_begin + evtRecEvent->totalCharged();

  EvtRecTrackIterator neutral_begin = track_begin + evtRecEvent->totalCharged();
  EvtRecTrackIterator neutral_end = track_begin + evtRecEvent->totalTracks();

  int runNo = eventHeader->runNumber();
  int event = eventHeader->eventNumber();
  Ncut[0]++;

  m_run = runNo;
  m_event = event;
  m_rec_run = runNo;
  m_rec_event = event;
  m_rec_nchrg = evtRecEvent->totalCharged();
  m_rec_nneu = evtRecEvent->totalNeutral();

  if (0 == (m_nEvt % m_nEvtDisp))
    std::cout << "Run num : " << runNo << ",  Event " << m_nEvt
              << ",   Event number online " << event << endl;
  m_nEvt++;

  // 2020.03.27 for topology
  ////////////////////////////////////////////////////////////////////////////////////wu
  ///begin
  //
  if (eventHeader->runNumber() < 0 && m_testMC == 1) {
    int m_numParticle = 0;
    double pow_tau = 3.1572315;
    double pow_r_2 = 0.0061699104;
    m_pi_e = 0.0;
    m_gam_e = 0.0;
    m_tau_e = 0.0;
    m_nu_e = 0.0;
    // MC information
    SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),
                                                     "/Event/MC/McParticleCol");
    if (!mcParticleCol) {
      std::cout << "Could not retrieve McParticelCol" << std::endl;
      return StatusCode::FAILURE;
    }
    Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();

    for (; iter_mc != mcParticleCol->end(); iter_mc++) {
      HepLorentzVector lvMC = (*iter_mc)->initialFourMomentum();
      //what do you mean by ma??
      ma_trackID[m_numParticle] = (*iter_mc)->particleProperty();
      ma_motherID[m_numParticle] = ((*iter_mc)->mother()).particleProperty();
      ma_trackIndex[m_numParticle] = (*iter_mc)->trackIndex();
      ma_motherIndex[m_numParticle] = ((*iter_mc)->mother()).trackIndex();
      ma_fromGenerator[m_numParticle] = (*iter_mc)->decayFromGenerator();
      ma_primaryParticle[m_numParticle] = (*iter_mc)->primaryParticle();
      ma_px[m_numParticle] = lvMC.px();
      ma_py[m_numParticle] = lvMC.py();
      ma_pz[m_numParticle] = lvMC.pz();
      ma_e[m_numParticle] = lvMC.e();

      if ((*iter_mc)->decayFromGenerator() &&
          ((*iter_mc)->mother()).particleProperty() == 15 &&
          lvMC.e() > m_pi_e) {
        m_pi_px = lvMC.px();
        m_pi_py = lvMC.py();
        m_pi_pz = lvMC.pz();
        m_pi_e = lvMC.e();
      }

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == 22 && lvMC.e() > m_gam_e) {
        m_gam_px = lvMC.px();
        m_gam_py = lvMC.py();
        m_gam_pz = lvMC.pz();
        m_gam_e = lvMC.e();
      }

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == 16 &&
          ((*iter_mc)->mother()).particleProperty() == 15 &&
          lvMC.e() > m_nu_e) {
        m_nu_px = lvMC.px();
        m_nu_py = lvMC.py();
        m_nu_pz = lvMC.pz();
        m_nu_e = lvMC.e();
      }

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == 15 && lvMC.e() > m_tau_e) {
        m_tau_px = lvMC.px();
        m_tau_py = lvMC.py();
        m_tau_pz = lvMC.pz();
        m_tau_e = lvMC.e();
      }
      m_numParticle += 1;
    }// end of loop of all particles
    
  
    HepLorentzVector p4Pi, p4Nu, p4Gam, p4Tau;
    setPxPyPzE(p4Pi, m_pi_px, m_pi_py, m_pi_pz, m_pi_e); 
    setPxPyPzE(p4Nu, m_nu_px, m_nu_py, m_nu_pz, m_nu_e); 
    setPxPyPzE(p4Gam, m_gam_px, m_gam_py, m_gam_pz, m_gam_e); 
    setPxPyPzE(p4Tau, m_tau_px, m_tau_py, m_tau_pz, m_tau_e); 

    m_x = 2.0 * p4Tau.dot(p4Gam) /pow_tau; 
    m_y = 2.0 * p4Tau.dot(p4Pi) /pow_tau; 
    m_z = m_x + m_y - 1.0; 
    m_t = (p4Pi+p4Gam).m2();  
    m_s = pow_tau * m_z;  //m_s is the same as m_t 

    m_within_kine_bounds = false; 
    if (m_x >= 0 && m_x <= (1.0 - pow_r_2) && m_y >= (1.0 - m_x + pow_r_2 / (1.0 - m_x)) && m_y <= (1 + pow_r_2)){
      m_within_kine_bounds= true; 
    }

    m_nTrack = m_numParticle;
  }// if this is MC and testMC is requested
  ///////////////////////////////////////////////////////////////////////////////////

  //***************************Access Truth information******************************
  // 2020.3.9, add the MC truth information
  if (eventHeader->runNumber() < 0) {
    mctruth();
    m_tuple1->write();

    m_tuple3->write();
  }

  if(m_writeGenOnly){
    return StatusCode::SUCCESS;
  }


  // reconstruct 2 good charged tracks, and don't limit the number of good
  // photons
  //***************************Select Good Charged
  //Tracks******************************
  // 2 good charged tracks
  Hep3Vector xorigin(0, 0, 0);
  VertexParameter m_privxpar;
  MyInitIP myInitIP;
  myInitIP.InitIP(xorigin, m_privxpar);

  Vint iGood;
  iGood.clear();
  Vint icharge;
  icharge.clear();
  int charge = 0;
  for (int i = 0; i < evtRecEvent->totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
    if (isGoodTrack(*itTrk, xorigin)) { // 2020.3.9, not use mdcKalTrack
      iGood.push_back(i);
      RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
      charge += mdcTrk->charge();
      icharge.push_back(mdcTrk->charge());
    }
  }
  int nGood = iGood.size();
  if (!(nGood == 2 && charge == 0))
    return StatusCode::SUCCESS;
  Ncut[1]++;

  // 2020.03.14, make sure the first track is positive, the second one is
  // negetive;

  if (icharge[0] == -1) {
    int temidx = iGood[0];
    iGood[0] = iGood[1];
    iGood[1] = temidx;
  }
  //***************************Select Good Photons******************************
  // nGood photons
  Vint iGam;
  iGam.clear();
  Vp4 pGam;
  pGam.clear();

  int igamma = 0;
  double esum = 0.;
  //The energy and id of gam with maximum energy 
  double tmp_e_gam_max = 0.; // 2020.04.13 modify from -100 to 0 !!
  int tmp_gam_max_id = -100;

  for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks(); i++) {
    if (i >= evtRecTrkCol->size()){
      break;
    }
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
    if (!isGoodShower(*itTrk)){
      continue;
    }
    RecEmcShower *aShower = NULL;
    aShower = (*itTrk)->emcShower();
    esum += aShower->energy();

    double tmp_e_gam = 0.;
    tmp_e_gam = aShower->energy();
    if (tmp_e_gam > tmp_e_gam_max) {
      tmp_e_gam_max = tmp_e_gam;
      tmp_gam_max_id = igamma;
    }

    MyInfoShower(aShower, m_shower_info, igamma);
    igamma++;
    iGam.push_back(i);
  }
  int nGam = iGam.size();
  m_nsig_shower = igamma;
  m_esum = esum;
  m_egam_max = tmp_e_gam_max;
  m_id_gam_max = tmp_gam_max_id;
  // 2020.03.16, save gamma track
  for (int i = 0; i < nGam; i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i];
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double eraw = emcTrk->energy();
    double phi = emcTrk->phi();
    double the = emcTrk->theta();
    HepLorentzVector ptrk;
    ptrk.setPx(eraw * sin(the) * cos(phi));
    ptrk.setPy(eraw * sin(the) * sin(phi));
    ptrk.setPz(eraw * cos(the));
    ptrk.setE(eraw);
    // ptrk = ptrk.boost(-0.011,0,0);// boost to cms
    pGam.push_back(ptrk);
  }

  int eventFlag = -1;
  // 2 or 3 good photons 
  if (nGam > 1 && nGam < 4) {
    eventFlag = 1;
  }
  // 4 or 5 good photons 
  if (nGam > 3 && nGam < 6) {
    eventFlag = 2;
  }

  if (igamma > 100){
    return StatusCode::SUCCESS;
  }
  Ncut[2]++;

  //***************************Save Information******************************
  // 2020.03.14-15, save information of charged tracks
  double theta_1;
  double theta_2;
  double phi_1;
  double phi_2;
  HepLorentzVector P4_1[5];
  HepLorentzVector P4_2[5];
  HepLorentzVector cms_P4_1[5];
  HepLorentzVector cms_P4_2[5];
  Hep3Vector P3_1;
  Hep3Vector P3_2;
  RecMdcKalTrack::PidType pidtype[5] = {
    RecMdcKalTrack::electron, RecMdcKalTrack::muon,  RecMdcKalTrack::pion,
    RecMdcKalTrack::kaon,     RecMdcKalTrack::proton
  };
  Vint itrk;
  for (int i = 0; i < nGood; i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
    if (!(*itTrk)->isMdcTrackValid()){
      continue;
    }
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
    if ((*itTrk)->isEmcShowerValid()) {
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      double pctrk = mdcTrk->p();
      double thetamdc = mdcTrk->theta();
      double eemc = emcTrk->energy();
      // The E/p      
      double evp = eemc / pctrk;
      itrk.push_back(iGood[i]);
      if (i < 1) {
        // Mdc information for the first track
        RecMdcKalTrack *mdckaltrk = (*itTrk)->mdcKalTrack();
        for (int j = 0; j < 5; j++) {
          RecMdcKalTrack::setPidType(pidtype[j]);
          HepLorentzVector ptrk;
          ptrk.setPx(mdckaltrk->px());
          ptrk.setPy(mdckaltrk->py());
          ptrk.setPz(mdckaltrk->pz());
          double p3 = ptrk.mag();
          ptrk.setE(sqrt(p3 * p3 + xmass[j] * xmass[j]));
          P4_1[j] = ptrk;
        }
        m_charge_1 = mdckaltrk->charge();
        theta_1 = mdcTrk->theta();
        phi_1 = mdcTrk->phi();

        m_pmdc_1 = mdcTrk->p();
        m_ptmdc_1 = mdcTrk->pxy();

        for (int j = 0; j < 5; j++) {
          m_pkal_1[j] = P4_1[j].rho();
          m_ptkal_1[j] = P4_1[j].perp();
        }

        m_theta_1 = thetamdc;
        m_p3_1[0] = mdcTrk->px();
        m_p3_1[1] = mdcTrk->py();
        m_p3_1[2] = mdcTrk->pz();
        P3_1.set(mdcTrk->px(), mdcTrk->py(), mdcTrk->pz());

        // PID information for the first track
        ParticleID *pid = ParticleID::instance();
        pid->init();
        pid->setMethod(pid->methodProbability());
        pid->setChiMinCut(4);
        pid->setRecTrack(*itTrk);
        pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() |
                       pid->useEmc());
        pid->identify(pid->onlyProton() | pid->onlyKaon() | pid->onlyPion() |
                      pid->onlyElectron() | pid->onlyMuon());
        pid->calculate();
        if (pid->IsPidInfoValid()) {
          m_probe_1 = pid->probElectron();
          m_probmu_1 = pid->probMuon();
          m_probpi_1 = pid->probPion();
          m_probk_1 = pid->probKaon();
          m_probp_1 = pid->probProton();
        } else {
          m_probe_1 = -1.0;
          m_probmu_1 = -1.0;
          m_probpi_1 = -1.0;
          m_probk_1 = -1.0;
          m_probp_1 = -1.0;
        }
        // TOF information for the first track
        if ((*itTrk)->isTofTrackValid()) {
          SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
          SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
          int nTofInfo = 0; 
          //@TODO: The following is wrong.  
          for (; iter_tof != tofTrkCol.end(); iter_tof++) {
            TofHitStatus *status = new TofHitStatus;
            status->setStatus((*iter_tof)->status());
            m_fgtof_1 = 1;
            m_counter_1 = status->is_counter();
            m_isbarrel_1 = status->is_barrel();
            m_layertof_1 = status->layer();
            m_iscluster_1 = status->is_cluster();
            m_tof_1 = (*iter_tof)->tof();
            m_texe_1 = (*iter_tof)->texpElectron();
            m_texmu_1 = (*iter_tof)->texpMuon();
            m_texpi_1 = (*iter_tof)->texpPion();
            m_texk_1 = (*iter_tof)->texpKaon();
            m_texp_1 = (*iter_tof)->texpProton();
            m_dte_1 = m_tof_1 - m_texe_1;
            m_dtmu_1 = m_tof_1 - m_texmu_1;
            m_dtpi_1 = m_tof_1 - m_texpi_1;
            m_dtk_1 = m_tof_1 - m_texk_1;
            m_dtp_1 = m_tof_1 - m_texp_1;
            nTofInfo ++ ;
          }
          //std::cout<<"There are " <<nTofInfo <<" Tof info"<< std::endl;
        }

        // EMC information for the first track
        m_evp_1 = evp;
        m_ene_1 = eemc;

        // dedx information for the first track
        if ((*itTrk)->isMdcDedxValid()) {
          RecMdcDedx *dedxTrk = (*itTrk)->mdcDedx();
          m_chi_e_1 = dedxTrk->chiE();
          m_chi_mu_1 = dedxTrk->chiMu();
          m_chi_pi_1 = dedxTrk->chiPi();
          m_chi_k_1 = dedxTrk->chiK();
          m_chi_p_1 = dedxTrk->chiP();
        } else {
          m_chi_e_1 = -100.;
          m_chi_mu_1 = -100.;
          m_chi_pi_1 = -100.;
          m_chi_k_1 = -100.;
          m_chi_p_1 = -100.;
        }

        // MUC information for the first track
        if ((*itTrk)->isMucTrackValid()) {
          RecMucTrack *mucTrk = (*itTrk)->mucTrack();
          m_maxhitsinlay_1 = mucTrk->maxHitsInLayer();
          m_numhits_1 = mucTrk->numHits();
          m_numlayers_1 = mucTrk->numLayers();
          m_depth_1 = mucTrk->depth();
          m_mucchi2_1 = mucTrk->chi2();
        } else {
          m_numhits_1 = -10;
          m_maxhitsinlay_1 = -10;
          m_numlayers_1 = -10;
          m_depth_1 = -10;
          m_mucchi2_1 = 300;
        }
      } // the end of the information of the first track
      else {
        // Mdc information for the second track
        RecMdcKalTrack *mdckaltrk = (*itTrk)->mdcKalTrack();
        for (int j = 0; j < 5; j++) {
          RecMdcKalTrack::setPidType(pidtype[j]);
          HepLorentzVector ptrk;
          ptrk.setPx(mdckaltrk->px());
          ptrk.setPy(mdckaltrk->py());
          ptrk.setPz(mdckaltrk->pz());
          double p3 = ptrk.mag();
          ptrk.setE(sqrt(p3 * p3 + xmass[j] * xmass[j]));
          P4_2[j] = ptrk;
        }

        m_charge_2 = mdckaltrk->charge();
        theta_2 = mdcTrk->theta();
        phi_2 = mdcTrk->phi();

        m_pmdc_2 = mdcTrk->p();
        m_ptmdc_2 = mdcTrk->pxy();
        for (int j = 0; j < 5; j++) {
          m_pkal_2[j] = P4_2[j].rho();
          m_ptkal_2[j] = P4_2[j].perp();
        }

        m_theta_2 = thetamdc;
        m_p3_2[0] = mdcTrk->px();
        m_p3_2[1] = mdcTrk->py();
        m_p3_2[2] = mdcTrk->pz();
        P3_2.set(mdcTrk->px(), mdcTrk->py(), mdcTrk->pz());

        // PID information for the second track
        ParticleID *pid = ParticleID::instance();
        pid->init();
        pid->setMethod(pid->methodProbability());
        pid->setChiMinCut(4);
        pid->setRecTrack(*itTrk);
        pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() |
                       pid->useEmc());
        pid->identify(pid->onlyProton() | pid->onlyKaon() | pid->onlyPion() |
                      pid->onlyElectron() | pid->onlyMuon());
        pid->calculate();
        if (pid->IsPidInfoValid()) {
          m_probe_2 = pid->probElectron();
          m_probmu_2 = pid->probMuon();
          m_probpi_2 = pid->probPion();
          m_probk_2 = pid->probKaon();
          m_probp_2 = pid->probProton();
        } else {
          m_probe_2 = -1.0;
          m_probmu_2 = -1.0;
          m_probpi_2 = -1.0;
          m_probk_2 = -1.0;
          m_probp_2 = -1.0;
        }

        // TOF information for the second track
        if ((*itTrk)->isTofTrackValid()) {
          SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
          SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
          //@TODO: The following is wrong.  
          for (; iter_tof != tofTrkCol.end(); iter_tof++) {
            TofHitStatus *status = new TofHitStatus;
            status->setStatus((*iter_tof)->status());
            m_fgtof_2 = 1;
            m_counter_2 = status->is_counter();
            m_isbarrel_2 = status->is_barrel();
            m_layertof_2 = status->layer();
            m_iscluster_2 = status->is_cluster();
            m_tof_2 = (*iter_tof)->tof();
            m_texe_2 = (*iter_tof)->texpElectron();
            m_texmu_2 = (*iter_tof)->texpMuon();
            m_texpi_2 = (*iter_tof)->texpPion();
            m_texk_2 = (*iter_tof)->texpKaon();
            m_texp_2 = (*iter_tof)->texpProton();
            m_dte_2 = m_tof_2 - m_texe_2;
            m_dtmu_2 = m_tof_2 - m_texmu_2;
            m_dtpi_2 = m_tof_2 - m_texpi_2;
            m_dtk_2 = m_tof_2 - m_texk_2;
            m_dtp_2 = m_tof_2 - m_texp_2;
          }
        }

        // EMC information for the second track
        m_evp_2 = evp;
        m_ene_2 = eemc;

        // dedx information for the second track
        if ((*itTrk)->isMdcDedxValid()) {
          RecMdcDedx *dedxTrk = (*itTrk)->mdcDedx();
          m_chi_e_2 = dedxTrk->chiE();
          m_chi_mu_2 = dedxTrk->chiMu();
          m_chi_pi_2 = dedxTrk->chiPi();
          m_chi_k_2 = dedxTrk->chiK();
          m_chi_p_2 = dedxTrk->chiP();
        } else {
          m_chi_e_2 = -100.;
          m_chi_mu_2 = -100.;
          m_chi_pi_2 = -100.;
          m_chi_k_2 = -100.;
          m_chi_p_2 = -100.;
        }

        // MUC information for the second track
        if ((*itTrk)->isMucTrackValid()) {
          RecMucTrack *mucTrk = (*itTrk)->mucTrack();
          m_maxhitsinlay_2 = mucTrk->maxHitsInLayer();
          m_numhits_2 = mucTrk->numHits();
          m_numlayers_2 = mucTrk->numLayers();
          m_depth_2 = mucTrk->depth();
          m_mucchi2_2 = mucTrk->chi2();
        } else {
          m_numhits_2 = -10;
          m_maxhitsinlay_2 = -10;
          m_numlayers_2 = -10;
          m_depth_2 = -10;
          m_mucchi2_2 = 300;
        }
      } // the end of the information of the second track
    }
  } // the end of the infotmation of charged tracks
  Ncut[3]++;

  // 2020.03.15, to make sure the different ID of two tracks
  int nchrgtrk = itrk.size();

  // 2020.04.24, check the cut - nTrack in EMC==2
  m_nTrk_EMC = nchrgtrk;
  if (nchrgtrk != 2){
    return StatusCode::SUCCESS;
  }
  Ncut[4]++;
  //How can this happen? 
  if (itrk[0] == itrk[1]){
    return StatusCode::SUCCESS;
  }
  Ncut[5]++;

  // 2020.03.15, to make sure the opposite charge of two tracks
  if (!(m_charge_1 * m_charge_2 < 0)){
    return StatusCode::SUCCESS;
  }
  Ncut[6]++;

  // 2020.03.15, calculate total momentum of all tracks and photons
  HepLorentzVector pkal_tot[5];
  HepLorentzVector pkal_Trks[5];
  for (int i = 0; i < 5; i++) {
    pkal_tot[i] = P4_1[i] + P4_2[i];
    pkal_Trks[i] = P4_1[i] + P4_2[i];
    for (int j = 0; j < nGam; j++) {
      pkal_tot[i] += pGam[j];
    }
  }

  for (int i = 0; i < 5; i++) {
    m_pkal_tot_px[i] = pkal_tot[i].px();
    m_pkal_tot_py[i] = pkal_tot[i].py();
    m_pkal_tot_pz[i] = pkal_tot[i].pz();
    m_pkal_tot_e[i] = pkal_tot[i].e();
    m_mTrks[i] = pkal_Trks[i].m();
    m_m2Trks[i] = pkal_Trks[i].m2();
  }

  m_ang_mdc_acol = P3_1.polarAngle(P3_2) * 180 / CLHEP::pi;
  m_ang_mdc_acop = P3_1.deltaPhi(P3_2) * 180 / CLHEP::pi;

  m_the_add = theta_1 + theta_2;
  m_phi_diff = phi_1 - phi_2;
  if (m_phi_diff < 0)
    m_phi_diff = m_phi_diff + twopi;
  m_PTEM = (P3_1 + P3_2).perp() / (m_ecms - (P3_1 + P3_2).r());

  HepLorentzVector ecms(0.011 * m_ecms, 0, 0, m_ecms);
  // boost to cms
  Hep3Vector P3_trk1[5];
  Hep3Vector P3_trk2[5];
  for (int i = 0; i < 5; i++) {
    /*		  P4_1[i].boost(-0.011,0,0);
                      P4_2[i].boost(-0.011,0,0);
                      cms_P4_1[i]=P4_1[i];
                      cms_P4_2[i]=P4_2[i];
                      P3_trk1[i] = P4_1[i].vect();
                      P3_trk2[i] = P4_2[i].vect();
                      //boost to lab system
                      P4_1[i].boost(0.011,0,0);
                      P4_2[i].boost(0.011,0,0);  */

    cms_P4_1[i] = P4_1[i];
    cms_P4_2[i] = P4_2[i];
    cms_P4_1[i].boost(-0.011, 0, 0);
    cms_P4_2[i].boost(-0.011, 0, 0);
    P3_trk1[i] = cms_P4_1[i].vect();
    P3_trk2[i] = cms_P4_2[i].vect();
  }
  HepLorentzVector p_all(0, 0, 0, m_ecms);
  HepLorentzVector P4_ee = cms_P4_1[0] + cms_P4_2[0];
  HepLorentzVector P4_emu = cms_P4_1[0] + cms_P4_2[1];
  HepLorentzVector P4_mue = cms_P4_1[1] + cms_P4_2[0];
  HepLorentzVector P4_epi = cms_P4_1[0] + cms_P4_2[2];
  HepLorentzVector P4_pie = cms_P4_1[2] + cms_P4_2[0];
  HepLorentzVector P4_ek = cms_P4_1[0] + cms_P4_2[3];
  HepLorentzVector P4_ke = cms_P4_1[3] + cms_P4_2[0];
  HepLorentzVector P4_mumu = cms_P4_1[1] + cms_P4_2[1];
  HepLorentzVector P4_mupi = cms_P4_1[1] + cms_P4_2[2];
  HepLorentzVector P4_pimu = cms_P4_1[2] + cms_P4_2[1];
  HepLorentzVector P4_pipi = cms_P4_1[2] + cms_P4_2[2];
  HepLorentzVector P4_muk = cms_P4_1[1] + cms_P4_2[3];
  HepLorentzVector P4_kmu = cms_P4_1[3] + cms_P4_2[1];
  HepLorentzVector P4_pik = cms_P4_1[2] + cms_P4_2[3];
  HepLorentzVector P4_kpi = cms_P4_1[3] + cms_P4_2[2];
  HepLorentzVector P4_kk = cms_P4_1[3] + cms_P4_2[3];
  m_miss_m2_ee = (m_ecms - P4_ee.e()) * (m_ecms - P4_ee.e()) -
                 ((p_all - P4_ee).rho()) * ((p_all - P4_ee).rho());
  m_miss_m2_emu = (m_ecms - P4_emu.e()) * (m_ecms - P4_emu.e()) -
                  ((p_all - P4_emu).rho()) * ((p_all - P4_emu).rho());
  m_miss_m2_mue = (m_ecms - P4_mue.e()) * (m_ecms - P4_mue.e()) -
                  ((p_all - P4_mue).rho()) * ((p_all - P4_mue).rho());
  m_miss_m2_epi = (m_ecms - P4_epi.e()) * (m_ecms - P4_epi.e()) -
                  ((p_all - P4_epi).rho()) * ((p_all - P4_epi).rho());
  m_miss_m2_pie = (m_ecms - P4_pie.e()) * (m_ecms - P4_pie.e()) -
                  ((p_all - P4_pie).rho()) * ((p_all - P4_pie).rho());
  m_miss_m2_ek = (m_ecms - P4_ek.e()) * (m_ecms - P4_ek.e()) -
                 ((p_all - P4_ek).rho()) * ((p_all - P4_ek).rho());
  m_miss_m2_ke = (m_ecms - P4_ke.e()) * (m_ecms - P4_ke.e()) -
                 ((p_all - P4_ke).rho()) * ((p_all - P4_ke).rho());
  m_miss_m2_mumu = (m_ecms - P4_mumu.e()) * (m_ecms - P4_mumu.e()) -
                   ((p_all - P4_mumu).rho()) * ((p_all - P4_mumu).rho());
  m_miss_m2_mupi = (m_ecms - P4_mupi.e()) * (m_ecms - P4_mupi.e()) -
                   ((p_all - P4_mupi).rho()) * ((p_all - P4_mupi).rho());
  m_miss_m2_pimu = (m_ecms - P4_pimu.e()) * (m_ecms - P4_pimu.e()) -
                   ((p_all - P4_pimu).rho()) * ((p_all - P4_pimu).rho());
  m_miss_m2_pipi = (m_ecms - P4_pipi.e()) * (m_ecms - P4_pipi.e()) -
                   ((p_all - P4_pipi).rho()) * ((p_all - P4_pipi).rho());
  m_miss_m2_muk = (m_ecms - P4_muk.e()) * (m_ecms - P4_muk.e()) -
                  ((p_all - P4_muk).rho()) * ((p_all - P4_muk).rho());
  m_miss_m2_kmu = (m_ecms - P4_kmu.e()) * (m_ecms - P4_kmu.e()) -
                  ((p_all - P4_kmu).rho()) * ((p_all - P4_kmu).rho());
  m_miss_m2_pik = (m_ecms - P4_pik.e()) * (m_ecms - P4_pik.e()) -
                  ((p_all - P4_pik).rho()) * ((p_all - P4_pik).rho());
  m_miss_m2_kpi = (m_ecms - P4_kpi.e()) * (m_ecms - P4_kpi.e()) -
                  ((p_all - P4_kpi).rho()) * ((p_all - P4_kpi).rho());
  m_miss_m2_kk = (m_ecms - P4_kk.e()) * (m_ecms - P4_kk.e()) -
                 ((p_all - P4_kk).rho()) * ((p_all - P4_kk).rho());

  m_ee_angle = P3_trk1[0].angle(P3_trk2[0]) * 180 / CLHEP::pi;
  m_emu_angle = P3_trk1[0].angle(P3_trk2[1]) * 180 / CLHEP::pi;
  m_mue_angle = P3_trk1[1].angle(P3_trk2[0]) * 180 / CLHEP::pi;
  m_epi_angle = P3_trk1[0].angle(P3_trk2[2]) * 180 / CLHEP::pi;
  m_pie_angle = P3_trk1[2].angle(P3_trk2[0]) * 180 / CLHEP::pi;
  m_ek_angle = P3_trk1[0].angle(P3_trk2[3]) * 180 / CLHEP::pi;
  m_ke_angle = P3_trk1[3].angle(P3_trk2[0]) * 180 / CLHEP::pi;
  m_mumu_angle = P3_trk1[1].angle(P3_trk2[1]) * 180 / CLHEP::pi;
  m_mupi_angle = P3_trk1[1].angle(P3_trk2[2]) * 180 / CLHEP::pi;
  m_pimu_angle = P3_trk1[2].angle(P3_trk2[1]) * 180 / CLHEP::pi;
  m_pipi_angle = P3_trk1[2].angle(P3_trk2[2]) * 180 / CLHEP::pi;
  m_muk_angle = P3_trk1[1].angle(P3_trk2[3]) * 180 / CLHEP::pi;
  m_kmu_angle = P3_trk1[3].angle(P3_trk2[1]) * 180 / CLHEP::pi;
  m_pik_angle = P3_trk1[2].angle(P3_trk2[3]) * 180 / CLHEP::pi;
  m_kpi_angle = P3_trk1[3].angle(P3_trk2[2]) * 180 / CLHEP::pi;
  m_kk_angle = P3_trk1[3].angle(P3_trk2[3]) * 180 / CLHEP::pi;

  m_ee_acop = P3_trk1[0].deltaPhi(P3_trk2[0]) * 180 / CLHEP::pi;
  m_emu_acop = P3_trk1[0].deltaPhi(P3_trk2[1]) * 180 / CLHEP::pi;
  m_mue_acop = P3_trk1[1].deltaPhi(P3_trk2[0]) * 180 / CLHEP::pi;
  m_epi_acop = P3_trk1[0].deltaPhi(P3_trk2[2]) * 180 / CLHEP::pi;
  m_pie_acop = P3_trk1[2].deltaPhi(P3_trk2[0]) * 180 / CLHEP::pi;
  m_ek_acop = P3_trk1[0].deltaPhi(P3_trk2[3]) * 180 / CLHEP::pi;
  m_ke_acop = P3_trk1[3].deltaPhi(P3_trk2[0]) * 180 / CLHEP::pi;
  m_mumu_acop = P3_trk1[1].deltaPhi(P3_trk2[1]) * 180 / CLHEP::pi;
  m_mupi_acop = P3_trk1[1].deltaPhi(P3_trk2[2]) * 180 / CLHEP::pi;
  m_pimu_acop = P3_trk1[2].deltaPhi(P3_trk2[1]) * 180 / CLHEP::pi;
  m_pipi_acop = P3_trk1[2].deltaPhi(P3_trk2[2]) * 180 / CLHEP::pi;
  m_muk_acop = P3_trk1[1].deltaPhi(P3_trk2[3]) * 180 / CLHEP::pi;
  m_kmu_acop = P3_trk1[3].deltaPhi(P3_trk2[1]) * 180 / CLHEP::pi;
  m_pik_acop = P3_trk1[2].deltaPhi(P3_trk2[3]) * 180 / CLHEP::pi;
  m_kpi_acop = P3_trk1[3].deltaPhi(P3_trk2[2]) * 180 / CLHEP::pi;
  m_kk_acop = P3_trk1[3].deltaPhi(P3_trk2[3]) * 180 / CLHEP::pi;

  m_ee_PTEM = (P3_trk1[0] + P3_trk2[0]).perp() /
              (m_ecms - (P3_trk1[0] + P3_trk2[0]).r());
  m_emu_PTEM = (P3_trk1[0] + P3_trk2[1]).perp() /
               (m_ecms - (P3_trk1[0] + P3_trk2[1]).r());
  m_mue_PTEM = (P3_trk1[1] + P3_trk2[0]).perp() /
               (m_ecms - (P3_trk1[1] + P3_trk2[0]).r());
  m_epi_PTEM = (P3_trk1[0] + P3_trk2[2]).perp() /
               (m_ecms - (P3_trk1[0] + P3_trk2[2]).r());
  m_pie_PTEM = (P3_trk1[2] + P3_trk2[0]).perp() /
               (m_ecms - (P3_trk1[2] + P3_trk2[0]).r());
  m_ek_PTEM = (P3_trk1[0] + P3_trk2[3]).perp() /
              (m_ecms - (P3_trk1[0] + P3_trk2[3]).r());
  m_ke_PTEM = (P3_trk1[3] + P3_trk2[0]).perp() /
              (m_ecms - (P3_trk1[3] + P3_trk2[0]).r());
  m_mumu_PTEM = (P3_trk1[1] + P3_trk2[1]).perp() /
                (m_ecms - (P3_trk1[1] + P3_trk2[1]).r());
  m_mupi_PTEM = (P3_trk1[1] + P3_trk2[2]).perp() /
                (m_ecms - (P3_trk1[1] + P3_trk2[2]).r());
  m_pimu_PTEM = (P3_trk1[2] + P3_trk2[1]).perp() /
                (m_ecms - (P3_trk1[2] + P3_trk2[1]).r());
  m_pipi_PTEM = (P3_trk1[2] + P3_trk2[2]).perp() /
                (m_ecms - (P3_trk1[2] + P3_trk2[2]).r());
  m_muk_PTEM = (P3_trk1[1] + P3_trk2[3]).perp() /
               (m_ecms - (P3_trk1[1] + P3_trk2[3]).r());
  m_kmu_PTEM = (P3_trk1[3] + P3_trk2[1]).perp() /
               (m_ecms - (P3_trk1[3] + P3_trk2[1]).r());
  m_pik_PTEM = (P3_trk1[2] + P3_trk2[3]).perp() /
               (m_ecms - (P3_trk1[2] + P3_trk2[3]).r());
  m_kpi_PTEM = (P3_trk1[3] + P3_trk2[2]).perp() /
               (m_ecms - (P3_trk1[3] + P3_trk2[2]).r());
  m_kk_PTEM = (P3_trk1[3] + P3_trk2[3]).perp() /
              (m_ecms - (P3_trk1[3] + P3_trk2[3]).r());

  // 2020.03.15, reconstruct pi0 and rho, in the variable, for example mrho_mn,
  // m is the index of pi0, n is the index of charged pi
  if (eventFlag == 1) { // pi0
    double delmpi0 = 999.;
    double mpi0 = 0.1349766;
    int g1 = -1;
    int g2 = -1;
    HepLorentzVector pTot;
    for (int i = 0; i < nGam - 1; i++) {
      for (int j = i + 1; j < nGam; j++) {
        HepLorentzVector p2g = pGam[i] + pGam[j];
        double m2g = p2g.m();
        double dltm2g = abs(m2g - mpi0);
        if (delmpi0 > dltm2g) {
          delmpi0 = dltm2g;
          g1 = i;
          g2 = j;
        }
      }
    }
    m_dlt_mpi0 = delmpi0;
    m_eg_11 = pGam[g1].e();
    m_eg_12 = pGam[g2].e();
    m_eg_21 = -1.0;
    m_eg_22 = -1.0;
    m_mpi0_1 = (pGam[g1] + pGam[g2]).m();
    m_mpi0_2 = -1.0;
    m_mrho_11 = (P4_1[2] + pGam[g1] + pGam[g2]).m();
    m_mrho_12 = (P4_2[2] + pGam[g1] + pGam[g2]).m();
    m_mrho_21 = -1.0;
    m_mrho_22 = -1.0;
    m_prho_11 = (P4_1[2] + pGam[g1] + pGam[g2]).rho();
    m_prho_12 = (P4_2[2] + pGam[g1] + pGam[g2]).rho();
    m_prho_21 = -10.0;
    m_prho_22 = -10.0;

    // 2020.03.15, calculate the helicity angle
    HepLorentzVector ppi0_1 = pGam[g1] + pGam[g2];
    Hep3Vector bv_pi0_rec = ppi0_1.boostVector();
    HepLorentzVector p4rho_11 = P4_1[2] + pGam[g1] + pGam[g2];
    HepLorentzVector p4rho_12 = P4_2[2] + pGam[g1] + pGam[g2];
    Hep3Vector bv_rho_rec_11 = p4rho_11.boostVector();
    Hep3Vector bv_rho_rec_12 = p4rho_12.boostVector();

    // 2020.03.15, calculate acop angle
    Hep3Vector P3_rho_11 = p4rho_11.vect();
    Hep3Vector P3_rho_12 = p4rho_12.vect();
    m_erho_11_acop = P3_trk2[0].deltaPhi(P3_rho_11) * 180 / CLHEP::pi;
    m_erho_12_acop = P3_trk1[0].deltaPhi(P3_rho_12) * 180 / CLHEP::pi;
    m_murho_11_acop = P3_trk2[1].deltaPhi(P3_rho_11) * 180 / CLHEP::pi;
    m_murho_12_acop = P3_trk1[1].deltaPhi(P3_rho_12) * 180 / CLHEP::pi;

    HepLorentzVector P4_erho_11 = cms_P4_2[0] + p4rho_11;
    HepLorentzVector P4_erho_12 = cms_P4_1[0] + p4rho_12;
    HepLorentzVector P4_murho_11 = cms_P4_2[1] + p4rho_11;
    HepLorentzVector P4_murho_12 = cms_P4_1[1] + p4rho_12;

    m_miss_m2_erho_11 =
        (m_ecms - P4_erho_11.e()) * (m_ecms - P4_erho_11.e()) -
        ((p_all - P4_erho_11).rho()) * ((p_all - P4_erho_11).rho());
    m_miss_m2_erho_12 =
        (m_ecms - P4_erho_12.e()) * (m_ecms - P4_erho_12.e()) -
        ((p_all - P4_erho_12).rho()) * ((p_all - P4_erho_12).rho());
    m_miss_m2_murho_11 =
        (m_ecms - P4_murho_11.e()) * (m_ecms - P4_murho_11.e()) -
        ((p_all - P4_murho_11).rho()) * ((p_all - P4_murho_11).rho());
    m_miss_m2_murho_12 =
        (m_ecms - P4_murho_12.e()) * (m_ecms - P4_murho_12.e()) -
        ((p_all - P4_murho_12).rho()) * ((p_all - P4_murho_12).rho());

    m_erho_11_angle = P3_trk2[0].angle(P3_rho_11) * 180 / CLHEP::pi;
    m_erho_12_angle = P3_trk1[0].angle(P3_rho_12) * 180 / CLHEP::pi;
    m_murho_11_angle = P3_trk2[1].angle(P3_rho_11) * 180 / CLHEP::pi;
    m_murho_12_angle = P3_trk1[1].angle(P3_rho_12) * 180 / CLHEP::pi;

    // 2020.03.15, calculate theta_1
    P4_1[2].boost(-bv_rho_rec_11);
    P4_2[2].boost(-bv_rho_rec_12);

    m_theta_pi_11 = P4_1[2].theta();
    m_theta_pi_12 = P4_2[2].theta();
    m_theta_pi_21 = -10.;
    m_theta_pi_22 = -10.;

    m_ppi_cms_11 = P4_1[2].rho();
    m_ppi_cms_12 = P4_2[2].rho();
    m_ppi_cms_21 = -10.;
    m_ppi_cms_22 = -10.;

    // 2020.03.15, calculate theta_2
    pGam[g1].boost(-bv_pi0_rec);
    pGam[g2].boost(-bv_pi0_rec);
    ppi0_1.boost(-bv_rho_rec_11);
    if (Ncut[6] % 2 == 0) {
      m_theta_m_11 = pGam[g1].vect().angle(ppi0_1.vect());
    } else {
      m_theta_m_11 = pGam[g2].vect().angle(ppi0_1.vect());
    }

    // 2020.03.15, boost pi0 to the second rho
    ppi0_1.boost(bv_rho_rec_11);
    ppi0_1.boost(-bv_rho_rec_12);

    if (Ncut[6] % 2 == 0) {
      m_theta_m_12 = pGam[g1].vect().angle(ppi0_1.vect());
    } else {
      m_theta_m_12 = pGam[g2].vect().angle(ppi0_1.vect());
    }
    m_theta_m_21 = -10.;
    m_theta_m_22 = -10.;
  } // end of reconstruct pi0

  if (eventFlag == 2) { // rho
    double delmpi0 = 999.;
    double mpi0 = 0.1349766;
    int g1 = -1;
    int g2 = -1;
    int g3 = -1;
    int g4 = -1;
    HepLorentzVector pTot;
    for (int i = 0; i < nGam - 3; i++) {
      for (int j = i + 1; j < nGam - 2; j++) {
        for (int k = j + 1; k < nGam - 1; k++) {
          for (int m = k + 1; m < nGam; m++) {
            HepLorentzVector p2g_1 = pGam[i] + pGam[j];
            HepLorentzVector p2g_2 = pGam[k] + pGam[m];
            double m2g_1 = p2g_1.m();
            double m2g_2 = p2g_2.m();
            double dltm2g = sqrt((m2g_1 - mpi0) * (m2g_1 - mpi0) +
                                 (m2g_2 - mpi0) * (m2g_2 - mpi0));
            if (delmpi0 > dltm2g) {
              delmpi0 = dltm2g;
              g1 = i;
              g2 = j;
              g3 = k;
              g4 = m;
            }
            p2g_1 = pGam[i] + pGam[k];
            p2g_2 = pGam[j] + pGam[m];
            m2g_1 = p2g_1.m();
            m2g_2 = p2g_2.m();
            dltm2g = sqrt((m2g_1 - mpi0) * (m2g_1 - mpi0) +
                          (m2g_2 - mpi0) * (m2g_2 - mpi0));
            if (delmpi0 > dltm2g) {
              delmpi0 = dltm2g;
              g1 = i;
              g2 = k;
              g3 = j;
              g4 = m;
            }
            p2g_1 = pGam[i] + pGam[m];
            p2g_2 = pGam[j] + pGam[k];
            m2g_1 = p2g_1.m();
            m2g_2 = p2g_2.m();
            dltm2g = sqrt((m2g_1 - mpi0) * (m2g_1 - mpi0) +
                          (m2g_2 - mpi0) * (m2g_2 - mpi0));
            if (delmpi0 > dltm2g) {
              delmpi0 = dltm2g;
              g1 = i;
              g2 = m;
              g3 = j;
              g4 = k;
            }
          }
        }
      }
    }
    m_dlt_mpi0 = delmpi0;
    m_eg_11 = pGam[g1].e();
    m_eg_12 = pGam[g2].e();
    m_eg_21 = pGam[g3].e();
    m_eg_22 = pGam[g4].e();
    m_mpi0_1 = (pGam[g1] + pGam[g2]).m();
    m_mpi0_2 = (pGam[g3] + pGam[g4]).m();
    m_mrho_11 = (P4_1[2] + pGam[g1] + pGam[g2]).m();
    m_mrho_12 = (P4_2[2] + pGam[g1] + pGam[g2]).m();
    m_mrho_21 = (P4_1[2] + pGam[g3] + pGam[g4]).m();
    m_mrho_22 = (P4_2[2] + pGam[g3] + pGam[g4]).m();
    m_prho_11 = (P4_1[2] + pGam[g1] + pGam[g2]).rho();
    m_prho_12 = (P4_2[2] + pGam[g1] + pGam[g2]).rho();
    m_prho_21 = (P4_1[2] + pGam[g3] + pGam[g4]).rho();
    m_prho_22 = (P4_2[2] + pGam[g3] + pGam[g4]).rho();

    // 2020.03.15, calculate the helicity angle
    // the first pi0
    HepLorentzVector ppi0_1 = pGam[g1] + pGam[g2];
    Hep3Vector bv_pi0_1_rec = ppi0_1.boostVector();
    HepLorentzVector p4rho_11 = P4_1[2] + pGam[g1] + pGam[g2];
    HepLorentzVector p4rho_12 = P4_2[2] + pGam[g1] + pGam[g2];
    Hep3Vector bv_rho_rec_11 = p4rho_11.boostVector();
    Hep3Vector bv_rho_rec_12 = p4rho_12.boostVector();

    // 2020.03.15, calculate acop angle
    Hep3Vector P3_rho_11 = p4rho_11.vect();
    Hep3Vector P3_rho_12 = p4rho_12.vect();
    // 2020.03.15, calculate theta_1
    P4_1[2].boost(-bv_rho_rec_11);
    P4_2[2].boost(-bv_rho_rec_12);

    m_theta_pi_11 = P4_1[2].theta();
    m_theta_pi_12 = P4_2[2].theta();

    m_ppi_cms_11 = P4_1[2].rho();
    m_ppi_cms_12 = P4_2[2].rho();

    // 2020.03.15, boost back to lab
    P4_1[2].boost(bv_rho_rec_11);
    P4_2[2].boost(bv_rho_rec_12);

    // 2020.03.15, the second pi0
    HepLorentzVector ppi0_2 = pGam[g3] + pGam[g4];
    Hep3Vector bv_pi0_2_rec = ppi0_2.boostVector();
    HepLorentzVector p4rho_21 = P4_1[2] + pGam[g3] + pGam[g4];
    HepLorentzVector p4rho_22 = P4_2[2] + pGam[g3] + pGam[g4];
    Hep3Vector bv_rho_rec_21 = p4rho_21.boostVector();
    Hep3Vector bv_rho_rec_22 = p4rho_22.boostVector();

    // 2020.03.15, calculate acop angle
    Hep3Vector P3_rho_21 = p4rho_21.vect();
    Hep3Vector P3_rho_22 = p4rho_22.vect();
    m_rhorho_11_acop = P3_rho_11.deltaPhi(P3_rho_22) * 180 / CLHEP::pi;
    m_rhorho_12_acop = P3_rho_12.deltaPhi(P3_rho_21) * 180 / CLHEP::pi;

    HepLorentzVector P4_rhorho_11 = p4rho_11 + p4rho_22;
    HepLorentzVector P4_rhorho_12 = p4rho_12 + p4rho_21;

    m_miss_m2_rhorho_11 =
        (m_ecms - P4_rhorho_11.e()) * (m_ecms - P4_rhorho_11.e()) -
        ((p_all - P4_rhorho_11).rho()) * ((p_all - P4_rhorho_11).rho());
    m_miss_m2_rhorho_12 =
        (m_ecms - P4_rhorho_12.e()) * (m_ecms - P4_rhorho_12.e()) -
        ((p_all - P4_rhorho_12).rho()) * ((p_all - P4_rhorho_12).rho());

    m_rhorho_11_angle = P3_rho_11.angle(P3_rho_22) * 180 / CLHEP::pi;
    m_rhorho_12_angle = P3_rho_12.angle(P3_rho_21) * 180 / CLHEP::pi;

    // 2020.03.15, calculate theta_1
    P4_1[2].boost(-bv_rho_rec_21);
    P4_2[2].boost(-bv_rho_rec_22);

    m_theta_pi_21 = P4_1[2].theta();
    m_theta_pi_22 = P4_2[2].theta();

    m_ppi_cms_21 = P4_1[2].rho();
    m_ppi_cms_22 = P4_2[2].rho();

    // 2020.03.15, calculate theta_2
    // the first pi0
    pGam[g1].boost(-bv_pi0_1_rec);
    pGam[g2].boost(-bv_pi0_1_rec);
    ppi0_1.boost(-bv_rho_rec_11);
    if (Ncut[6] % 2 == 0) {
      m_theta_m_11 = pGam[g1].vect().angle(ppi0_1.vect());
    } else {
      m_theta_m_11 = pGam[g2].vect().angle(ppi0_1.vect());
    }

    // 2020.03.15, boost pi0 to the second rho
    ppi0_1.boost(bv_rho_rec_11);
    ppi0_1.boost(-bv_rho_rec_12);

    if (Ncut[6] % 2 == 0) {
      m_theta_m_12 = pGam[g1].vect().angle(ppi0_1.vect());
    } else {
      m_theta_m_12 = pGam[g2].vect().angle(ppi0_1.vect());
    }

    // the second pi0
    pGam[g3].boost(-bv_pi0_2_rec);
    pGam[g4].boost(-bv_pi0_2_rec);
    ppi0_2.boost(-bv_rho_rec_21);
    if (Ncut[6] % 2 == 0) {
      m_theta_m_21 = pGam[g3].vect().angle(ppi0_2.vect());
    } else {
      m_theta_m_21 = pGam[g4].vect().angle(ppi0_2.vect());
    }

    // 2020.03.15, boost pi0 to the second rho
    ppi0_2.boost(bv_rho_rec_21);
    ppi0_2.boost(-bv_rho_rec_22);

    if (Ncut[6] % 2 == 0) {
      m_theta_m_22 = pGam[g3].vect().angle(ppi0_2.vect());
    } else {
      m_theta_m_22 = pGam[g4].vect().angle(ppi0_2.vect());
    }
  }
  if (eventFlag != 1 && eventFlag != 2) {
    m_eg_11 = -10.;
    m_eg_12 = -10.;
    m_eg_21 = -10.;
    m_eg_22 = -10.;
    m_mpi0_1 = -10.;
    m_mpi0_2 = -10.;
    m_mrho_11 = -10.;
    m_mrho_12 = -10.;
    m_mrho_21 = -10.;
    m_mrho_22 = -10.;
    m_prho_11 = -10.;
    m_prho_12 = -10.;
    m_prho_21 = -10.;
    m_prho_22 = -10.;
    m_theta_pi_11 = -10.;
    m_theta_pi_12 = -10.;
    m_theta_pi_21 = -10.;
    m_theta_pi_22 = -10.;
    m_theta_m_11 = -10.;
    m_theta_m_12 = -10.;
    m_theta_m_21 = -10.;
    m_theta_m_22 = -10.;
  }

  /*
  //2020.3.2, save information of charged tracks
  HepLorentzVector p4elp(0., 0., 0., 0.);
  HepLorentzVector p4mum(0., 0., 0., 0.);
  RecMdcKalTrack* aelpTrk       = NULL;
  RecMdcKalTrack* amumTrk       = NULL;
  int elp_id = -10;
  int mum_id = -10;
for(int k=0; k < iGood.size(); k++){
          EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[k];
          RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
          if(!(*itTrk)->isMdcTrackValid()) continue;
ParticleID *pid = ParticleID::instance();
pid->init();
pid->setMethod(pid->methodProbability());
pid->setChiMinCut(4);
pid->setRecTrack(*itTrk);
pid->usePidSys(pid->useDedx() | pid->useTof1() |  pid->useTof2() | pid->useEmc()
);
pid->identify(pid->onlyProton() | pid->onlyKaon()  | pid->onlyPion() |
pid->onlyElectron() | pid->onlyMuon());
pid->calculate();
if(!(pid->IsPidInfoValid()))     continue;

          double prob_electron      = -1.;
          double prob_muon          = -1.;
          double prob_pion          = -1.;
          double prob_kaon          = -1.;
          double prob_proton        = -1.;
          double prob_ratio_muon    = -1.;
          double prob_ratio_electron= -1.;
          prob_electron             = pid->probElectron();
          prob_muon                 = pid->probMuon();
          prob_pion                 = pid->probPion();
          prob_kaon                 = pid->probKaon();
          prob_proton               = pid->probProton();

prob_ratio_muon    = prob_muon / (prob_muon + prob_electron + prob_kaon);
prob_ratio_electron= prob_electron / (prob_pion + prob_electron + prob_kaon);

          if(prob_electron > 0.001 &&prob_electron > prob_pion && prob_electron
> prob_kaon){
                  TrackInfo track1Info;
                  track1Info.settrack(0, *itTrk);//electron information
                  MyInfoCharge(track1Info, m_elp_info);
m_elp_prob[0] = prob_electron ;
m_elp_prob[1] = prob_muon ;
m_elp_prob[2] = prob_pion ;
m_elp_prob[3] = prob_kaon ;
m_elp_prob[4] = prob_proton ;
m_elp_prob[5] = prob_ratio_electron ;
                  p4elp         = track1Info.p4();
      elp_id        = k;
}else if(prob_muon>0.001 &&prob_muon > prob_electron && prob_muon > prob_kaon){
                  TrackInfo track2Info;
                  track2Info.settrack(1, *itTrk);//muon information
                  MyInfoCharge(track2Info, m_mum_info);
m_mum_prob[0] = prob_electron ;
m_mum_prob[1] = prob_muon ;
m_mum_prob[2] = prob_pion ;
m_mum_prob[3] = prob_kaon ;
m_mum_prob[4] = prob_proton ;
m_mum_prob[5] = prob_ratio_muon ;
p4mum         =  track2Info.p4();
      mum_id        = k;
          }
  }

Ncut[3]++;
//2020.3.9, to make sure the opposite charge of electron and muon;
  if(!(m_elp_info[25] * m_mum_info[25] <0 )) return StatusCode::SUCCESS;
Ncut[4]++;
  if(!((elp_id >-1 && mum_id>-1) &&(elp_id + mum_id) == 1))  return
StatusCode::SUCCESS;
Ncut[5]++;
*/

  // 2020.3.9, save pi0 -> gamma e+ e-, eta->gamma e+ e-
  /*  CDPhotonList      photonList(neutral_begin, neutral_end, photonSelector);
          checkPi0Dalitz(photonList, p4elp, p4mum, metapdg, m_eta_gamee);
          checkPi0Dalitz(photonList, p4elp, p4mum, mpi0pdg, m_pi0_gamee);

    //2020.3.9, save gamma -> e+ e-
          RecMdcKalTrack *elpTrk =
    (*(evtRecTrkCol->begin()+iGood[elp_id]))->mdcKalTrack();
          RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
          RecMdcKalTrack *mumTrk =
    (*(evtRecTrkCol->begin()+iGood[mum_id]))->mdcKalTrack();
          //RecMdcKalTrack::setPidType(RecMdcKalTrack::muon);
          RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);//to see if there
    is electron misidentified as muon
    checkGammConversion(xorigin, elpTrk, mumTrk, m_gamcon);

    //2020.3.9, correction for electron due to bremsstrahlung, add the energy of
    photon from electron back to electron;
          HepLorentzVector p4elp_corfsr(0., 0., 0., 0.); p4elp_corfsr  = p4elp;
          std::vector<HepLorentzVector> p4FSRCor;
          p4FSRCor.clear();
          if (m_appFSRCorrection ) {
                  correctLeptonMomentumWithFSRPhoton(p4elp_corfsr, p4FSRCor);
                  p4elp_corfsr      += p4FSRCor[0];
                  setV4(m_gam_fromelp,    p4FSRCor[0]);// the total photons from
    electron
                  setV4(m_elp_corfsr,     p4elp_corfsr);// the electron after
    adding the photons from elctron
    }

          HepLorentzVector p4tot  (0., 0., 0., 0.);
          HepLorentzVector p4diff1(0., 0., 0., 0.);
          HepLorentzVector p4diff2(0., 0., 0., 0.);
          p4tot   = p4elp + p4mum;
          p4diff1 = p4mum - p4elp;
          p4diff2 = p4mum - p4elp_corfsr;
          m_info_elp_mum[0] = p4elp.vect().angle(p4mum.vect());
          m_info_elp_mum[1] = p4diff1.e() - p4diff1.rho();
          m_info_elp_mum[2] = p4diff2.e() - p4diff2.rho();
          setV4(m_info_emu,    p4tot,   0);// the total photons from electron
          setV4(m_info_emu,    p4diff1, 1);// the total photons from electron
          setV4(m_info_emu,    p4diff2, 2);// the total photons from electron
  */

  m_tuple2->write();
  return StatusCode::SUCCESS;
}

// Thanks for Zhou Xingyu's help 2020.03.28
bool tautauAlg::mctruth() {

  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),
                                                   "/Event/MC/McParticleCol");
  if (!mcParticleCol) {
    std::cout << "Could not retrieve McParticelCol" << std::endl;
    return (StatusCode::FAILURE);
  }

  // -------------------------------DDbar-----------------------------------
  Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
  const int incpdgid0 = 9030443; // 9030443 is the PDG code of Y(4260).
  m_idxmc = 0;
  m_rec_idxmc = 0;
  int ngamma = 0;
  bool incPdcy(false);
  int rootIndex(-1);
  for (; iter_mc != mcParticleCol->end(); iter_mc++) {
    // std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdg\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
    // if((*iter_mc)->primaryParticle()) continue;
    if (!(*iter_mc)->decayFromGenerator())
      continue;
    if (incPdcy == false && (*iter_mc)->particleProperty() == -22) {
      m_pdgid[m_idxmc] = 22222222; // 22222222 is the default PDG code of ISR
                                   // photons in TopoAna.
      m_rec_pdgid[m_idxmc] = 22222222; // 22222222 is the default PDG code of
                                       // ISR photons in TopoAna.
      m_motheridx[m_idxmc] = -1;
      m_rec_motheridx[m_idxmc] = -1;
      m_idxmc++;
      m_rec_idxmc++;
      ngamma++;
      continue;
    } else if ((*iter_mc)->particleProperty() == incpdgid0) {
      incPdcy = true;
      rootIndex = (*iter_mc)->trackIndex();
    } else if (incPdcy == true && (*iter_mc)->particleProperty() == 11 &&
               ((*iter_mc)->mother()).trackIndex() ==
                   (*iter_mc)->trackIndex()) {
      continue; // This statement is used to remove the redundant electron
                // following Y(4260).
    }
    if (!incPdcy)
      continue;
    m_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
    m_rec_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
    if ((*iter_mc)->particleProperty() == incpdgid0) {
      m_motheridx[m_idxmc] = -1;
      m_rec_motheridx[m_idxmc] = -1;
    } else if (((*iter_mc)->mother()).particleProperty() == incpdgid0) {
      m_motheridx[m_idxmc] =
          ((*iter_mc)->mother()).trackIndex() - rootIndex + ngamma;
      m_rec_motheridx[m_idxmc] =
          ((*iter_mc)->mother()).trackIndex() - rootIndex + ngamma;
    } else {
      m_motheridx[m_idxmc] =
          ((*iter_mc)->mother()).trackIndex() - rootIndex - 1 + ngamma;
      m_rec_motheridx[m_idxmc] =
          ((*iter_mc)->mother()).trackIndex() - rootIndex - 1 + ngamma;
    }
    m_idxmc++;
    m_rec_idxmc++;
  }

  // ------------------------------(hadrons,
  // BesTwogamma)----------------------------
  if (incPdcy == false) {
    iter_mc = mcParticleCol->begin();
    const int incpdgid1 = 91; // 91 is the PDG code of cluster
    const int incpdgid2 = 92; // 92 is the PDG code of string
    m_idxmc = 0;
    m_rec_idxmc = 0;
    ngamma = 0;
    incPdcy = false;
    rootIndex = -1;
    bool findTFRE = false; // TFRE is short for is The First Redundant Electron.
    for (; iter_mc != mcParticleCol->end(); iter_mc++) {
      // std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdg\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
      // if((*iter_mc)->primaryParticle()) continue;
      if (!(*iter_mc)->decayFromGenerator())
        continue;
      if (incPdcy == false && (*iter_mc)->particleProperty() == 22) {
        m_pdgid[m_idxmc] = 22222222;     // 22222222 is the PDG code of cluster.
        m_rec_pdgid[m_idxmc] = 22222222; // 22222222 is the PDG code of cluster.
        m_motheridx[m_idxmc] = -1;
        m_rec_motheridx[m_idxmc] = -1;
        m_idxmc++;
        m_rec_idxmc++;
        ngamma++;
        continue;
      } else if ((*iter_mc)->particleProperty() == incpdgid1 ||
                 (*iter_mc)->particleProperty() == incpdgid2) {
        incPdcy = true;
        rootIndex = (*iter_mc)->trackIndex();
      } else if (incPdcy == true && (*iter_mc)->particleProperty() == 11 &&
                 ((*iter_mc)->mother()).trackIndex() ==
                     (*iter_mc)->trackIndex()) {
        if (findTFRE == false) {
          m_pdgid[m_idxmc] = 11;
          m_rec_pdgid[m_idxmc] = 11;
          findTFRE = true;
        } else {
          m_pdgid[m_idxmc] = -11;
          m_rec_pdgid[m_idxmc] = -11;
        }
        m_motheridx[m_idxmc] = -1;
        m_rec_motheridx[m_idxmc] = -1;
        m_idxmc++;
        m_rec_idxmc++;
        continue;
      }
      if (!incPdcy || (*iter_mc)->particleProperty() == incpdgid1 ||
          (*iter_mc)->particleProperty() == incpdgid2)
        continue;
      m_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
      m_rec_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
      if (((*iter_mc)->mother()).particleProperty() == incpdgid1 ||
          ((*iter_mc)->mother()).particleProperty() == incpdgid2) {
        m_motheridx[m_idxmc] = -1;
        m_rec_motheridx[m_idxmc] = -1;
      } else {
        m_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - rootIndex - 1 + ngamma;
        m_rec_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - rootIndex - 1 + ngamma;
      }
      m_idxmc++;
      m_rec_idxmc++;
    }
  }

  // ---------------QED (bhabha, gamgam, mumu, tautau, qqbar,
  // qqbar2)------------------
  if (incPdcy == false) {
    iter_mc = mcParticleCol->begin();
    int pdgid = (*iter_mc)->particleProperty();
    m_idxmc = 0;
    m_rec_idxmc = 0;
    int idxOfTRE = -1; // TRE is short for The Redundant Electron. Statements
                       // related to the variable is used to remove the
                       // redundant electron in the tautau sample.
    int nOfTRQG = 0;   // TRQG is short for The Redundant quarks and gluons.
                     // Statements related to the variable is used to remove the
                     // d, dbar, u, ubar, s, sbar, c, cbar, and g in the qqbar
                     // and qqbar2 samples.
    for (; iter_mc != mcParticleCol->end(); iter_mc++) {
      // if(!(*iter_mc)->decayFromGenerator()) continue;
      // std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdg\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
      if (pdgid == 23 && (*iter_mc)->particleProperty() == 11 &&
          ((*iter_mc)->mother()).trackIndex() == m_idxmc) {
        idxOfTRE = m_idxmc;
        continue;
      }
      if (abs((*iter_mc)->particleProperty()) == 1 ||
          abs((*iter_mc)->particleProperty()) == 2 ||
          abs((*iter_mc)->particleProperty()) == 3 ||
          abs((*iter_mc)->particleProperty()) == 4 ||
          abs((*iter_mc)->particleProperty()) == 21) {
        nOfTRQG++;
        continue;
      }
      m_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
      m_rec_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
      if (idxOfTRE != -1 && ((*iter_mc)->mother()).trackIndex() > idxOfTRE) {
        m_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - 1 - nOfTRQG;
        m_rec_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - 1 - nOfTRQG;
      } else {
        m_motheridx[m_idxmc] = ((*iter_mc)->mother()).trackIndex() - nOfTRQG;
        m_rec_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - nOfTRQG;
      }
      m_idxmc++;
      m_rec_idxmc++;
    }
  }

  if (!mcParticleCol) {
    return false;
  } else {
    Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
    for (; iter_mc != mcParticleCol->end(); iter_mc++) {
      if (!(*iter_mc)->decayFromGenerator())
        continue;
      HepLorentzVector p_mc = (*iter_mc)->initialFourMomentum();

      int self_id = (*iter_mc)->particleProperty();
      int mother_id = ((*iter_mc)->mother()).particleProperty();
      int grandmother_id = (((*iter_mc)->mother()).mother()).particleProperty();

      if (self_id == 11 && mother_id == 15) { // tau^- -> e^- nu anti-nu
        setV4(m_info_tru, p_mc, 1);
      }
      if (self_id == -11 && mother_id == -15) { // tau^+ -> e^+ nu anti-nu
        setV4(m_info_tru, p_mc, 2);
      }
      if (self_id == -13 && mother_id == -15) { // tau^+ -> mu^+ nu anti-nu
        setV4(m_info_tru, p_mc, 3);
      }
      if (self_id == 13 && mother_id == 15) { // tau^- -> mu^- nu anti-nu
        setV4(m_info_tru, p_mc, 4);
      }
      if (self_id == 15) { // tau^-
        setV4(m_info_tru, p_mc, 5);
      }
      if (self_id == -15) { // tau^+
        setV4(m_info_tru, p_mc, 6);
      }
    }
  }
  return true;
}

bool tautauAlg::correctLeptonMomentumWithFSRPhoton(
    HepLorentzVector p4lepp, std::vector<HepLorentzVector> &p4FSRCor) {
  p4FSRCor.clear();

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(),
                                        "/Event/EvtRec/EvtRecEvent");
  if (!evtRecEvent) {
    std::cout << "Could not find EvtRecEvent" << std::endl;
    return false;
  }
  SmartDataPtr<EvtRecTrackCol> evtRecTrackCol(eventSvc(),
                                              "/Event/EvtRec/EvtRecTrackCol");
  if (!evtRecTrackCol) {
    std::cout << "Could not find EvtRecTrackCol" << std::endl;
    return false;
  }

  HepLorentzVector p4FSRLepp(0., 0., 0., 0.);
  for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks();
       i++) {
    EvtRecTrackIterator itTrk = evtRecTrackCol->begin() + i;
    if (!(*itTrk)->isEmcShowerValid())
      continue;
    if (!isGoodShower(*itTrk))
      continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double eraw = emcTrk->energy();
    double phi = emcTrk->phi();
    double the = emcTrk->theta();
    HepLorentzVector p4shower(eraw * sin(the) * cos(phi),
                              eraw * sin(the) * sin(phi), eraw * cos(the),
                              eraw);
    double dthelepp = 200.;
    dthelepp = p4lepp.vect().theta(p4shower.vect());
    dthelepp =
        fabs(fmod(dthelepp + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi);
    if (dthelepp < CLHEP::pi * 5. / 180.) {
      p4FSRLepp += p4shower;
      m_eleppFSRArray[m_nleppFSR] = p4shower.e();
      m_dtheleppFSRArray[m_nleppFSR] = dthelepp;
      m_nleppFSR++;
    }
  }
  p4FSRCor.push_back(p4FSRLepp);
  return true;
}

//***************************Output******************************
StatusCode tautauAlg::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
  cout << "Ntot=             " << Ncut[0] << endl;
  cout << "2_goodTrack=      " << Ncut[1] << endl;
  cout << "nGoodgam=         " << Ncut[2] << endl;
  cout << "nTrack0=          " << Ncut[3] << endl;
  cout << "2_track_emc=      " << Ncut[4] << endl;
  cout << "differet_id=      " << Ncut[5] << endl;
  cout << "opposite_charge=  " << Ncut[6] << endl;
  return StatusCode::SUCCESS;
}
