use crate::p01::ProcEngine;
use bincode::{Decode, Encode};
//use phf::phf_map;
use phf_macros::phf_map;
use regex::Regex;
use std::collections::HashMap;
use std::collections::HashSet;
use std::error::Error;

#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct Pea {
    pub aream: HashMap<String, PeaArea>,
}
#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct PeaArea {
    pub arid: String,
    pub provm: HashMap<String, PeaProv>,
}
#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct PeaProv {
    pub pvnm: String,
    pub gppv: f32,
    pub evpc: f32,
    pub subm: HashMap<String, PeaSub>,
}
use crate::p03::SubLoadProfRepr;

#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct PeaSub {
    pub sbid: String,
    pub feedm: HashMap<String, PeaFeed>,
    pub name: String,
    pub enam: String,
    pub area: String,
    pub arid: String,
    pub volt: String,
    pub cate: String,
    pub egat: String,
    pub state: String,
    pub conf: String,
    pub trax: String,
    pub mvax: String,
    pub feed: String,
    pub feno: usize,
    pub feeders: Vec<String>,
    pub trxn: usize,
    pub mvxn: i32,
    pub prov: String,
    pub sbtp: String,
    pub n1d_s: u64,
    pub n1d_f: u64,
    pub lp_rep_23: SubLoadProfRepr,
    pub lp_rep_24: SubLoadProfRepr,
    pub vspps: Vec<usize>,
    pub spps: Vec<usize>,
    pub repls: Vec<usize>,
    pub aojv: Vec<AojObj>,
}

#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct EvDistCalc {
    pub id: String,
    pub ev_no: f32,
    pub ev_pc: f32,
    pub ev_ds: f32,
}

#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct PeaFeed {
    pub fdid: String,
    pub tranm: HashMap<u64, PeaTrans>,
}
#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct PeaTrans {
    pub mets: Vec<PeaMeter>,
    pub trid: String,
    pub pea: String,
    pub n1d: u64,
    pub n1d_f: u64,
    pub ix: usize,
    pub lix: usize,
    pub own: String,
    pub mts: Vec<usize>,
    pub aojs: Vec<usize>,
    pub amps: Vec<usize>,
    pub muns: Vec<usize>,
    pub zons: Vec<usize>,
    pub sols: Vec<usize>,
    pub vols: Vec<usize>,
    pub vopw: f32,
    pub vose: f32,
    pub kw: f32,

    pub tr_tag: Option<String>,
    pub tr_fid: Option<String>,
    pub tr_lt: Option<f32>,
    pub tr_ln: Option<f32>,
    pub tr_cd: Option<f32>,
    pub tr_aoj: Option<String>,
    pub tr_pea: Option<String>,
    pub tr_kva: Option<f32>,
    pub tr_own: Option<String>,
    pub tr_loc: Option<String>,
    pub tr_n1d: Option<u64>,
    //pub ar: String,
    //pub ly: String,
    //pub ix: usize,
}
#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct PeaMeter {
    pub mt_ins: Option<String>,
    pub mt_pea: Option<String>,
    pub mt_tag: Option<String>,
    pub mt_phs: Option<String>,
    pub mt_x: Option<f32>,
    pub mt_y: Option<f32>,
    pub mt_lt: Option<f32>,
    pub mt_ln: Option<f32>,
    pub mt_aoj: Option<String>,
    pub tr_tag: Option<String>,
    pub tr_fid: Option<String>,
    pub tr_lt: Option<f32>,
    pub tr_ln: Option<f32>,
    pub tr_cd: Option<f32>,
    pub tr_aoj: Option<String>,
    pub tr_pea: Option<String>,
    pub tr_kva: Option<f32>,
    pub tr_own: Option<String>,
    pub tr_loc: Option<String>,
    pub tr_n1d: Option<u64>,
    pub mt_n1d: Option<u64>,
    pub ar: String,
    pub ly: String,
    pub ix: usize,
    //pub bills: Vec<PeaBill>,
    pub trsg: String,
    pub pea: String,
    pub ca: String,
    pub inst: String,
    pub rate: String,
    pub volt: String,
    pub mru: String,
    pub mat: String,
    pub main: String,
    pub kwh15: f32,
    pub kwh18: f32,
    pub amt19: f32,
    pub idx: usize,
    pub meth: i32,
    pub met_type: MeterAccType,
}

#[derive(Encode, Decode, Debug, Clone, Default)]
pub enum MeterAccType {
    #[default]
    Small,
    Large,
}

#[derive(Encode, Decode, Debug, Clone, Default)]
pub enum GridLevel {
    #[default]
    DisTrans,
    Feeder,
    Sub,
    Area,
}

#[derive(Encode, Decode, Debug, Clone, Default)]
pub enum SumType {
    #[default]
    Sum,
    Max,
    Min,
}

use strum::IntoEnumIterator;
use strum_macros::EnumIter;

#[derive(Encode, Decode, Debug, Clone, Default, EnumIter)]
pub enum VarType {
    #[default]
    None,
    NewCarRegVp01,
    GppVp02,
    MaxPosPowSubVs01,
    MaxNegPowSubVs02,
    VsppMvVs03,
    SppHvVs04,
    BigLotMvVs05,
    BigLotHvVs06,
    SubPowCapVs07,
    MaxPosPowFeederVf01,
    MaxNegPowFeederVf02,
    MaxPosDiffFeederVf03,
    MaxNegDiffFeederVf04,
    NoMeterTransVt01,
    SmallSellTrVt02,
    ChgStnCapTrVt03,
    ChgStnSellTrVt04,
    PwCapTriVt05,
    ZoneTrVt06,
    PopTrVt07,
    UnbalPowTrVt08,
    PkPowTrVt09,
    LargeSellTrVt10,
    AllNoMeterTrVt11,
    NoTrVt12,
    HmChgEvTrVc01,
    LvPowSatTrVc02,
    CntLvPowSatTrVc03,
    ChgStnCapVc04,
    ChgStnSellVc05,
    MvPowSatTrVc06,
    PowSolarVc07,
    MvVsppVc08,
    HvSppVc09,
    SmallSellVc10,
    LargeSellVc11,
    UnbalPowVc12,
    CntUnbalPowVc13,
    Uc1ValVc14,
    Uc2ValVc15,
    Uc3ValVc16,
    Uc1RankVc17,
    Uc2RankVc18,
    Uc3RankVc19,

    NoHmChgEvTr,
    PowHmChgEvTr,

    PkSelPowPhsAKw,
    PkSelPowPhsBKw,
    PkSelPowPhsCKw,
    PkSelPowPhsAvg,
    PkSelPowPhsMax,
    UnbalPowRate,
    TransLossKw,
    UnbalPowLossKw,
    CntTrUnbalLoss,
    CntTrSatLoss,
    TakeNote,
    /// How likely the province to have EV car
    EvCarLikely,
    /// How likely the province to be select
    SelectLikely,
    SubSolarPeekMw,
    SubSolarEnergy,
    SolarEnergy,
}

#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct AssVar {
    pub v: f32,
    pub l: GridLevel,
    pub t: VarType,
    pub s: SumType,
}

impl AssVar {
    pub fn val(v: f32) -> AssVar {
        AssVar {
            t: VarType::None,
            s: SumType::Sum,
            v,
            ..Default::default()
        }
    }
    pub fn new(t: VarType, s: SumType) -> AssVar {
        AssVar {
            t,
            s,
            ..Default::default()
        }
    }
}
#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct PeaAssVar {
    pub arid: String,
    pub pvid: String,
    pub sbid: String,
    pub fdid: String,
    pub n1d: u64,
    pub own: String,
    pub peano: String,
    pub aoj: String,
    pub set: u32,
    pub v: Vec<AssVar>,
    pub res: f32,
}

pub fn z2o(v: f32) -> f32 {
    if v == 0f32 {
        1f32
    } else {
        v
    }
}

impl PeaAssVar {
    pub fn div(&mut self, o: f32) {
        if o == 0f32 {
            return;
        }
        for v in self.v.iter_mut() {
            v.v /= o;
        }
    }

    pub fn nor(&mut self, o: &PeaAssVar) {
        for (v, o) in self.v.iter_mut().zip(o.v.iter()) {
            v.v /= z2o(o.v);
        }
    }

    pub fn copy(&mut self, o: &PeaAssVar, t: VarType) {
        self.v[t.clone() as usize].v = o.v[t.clone() as usize].v;
    }

    pub fn add(&mut self, o: &PeaAssVar) {
        for (v, o) in self.v.iter_mut().zip(o.v.iter()) {
            match v.s {
                SumType::Sum => v.v += o.v,
                SumType::Max => v.v = v.v.max(o.v),
                SumType::Min => v.v = v.v.min(o.v),
            }
        }
    }
    pub fn max(&mut self, o: &PeaAssVar) {
        if self.set == 0 {
            for (v, o) in self.v.iter_mut().zip(o.v.iter()) {
                v.v = o.v;
            }
            return;
        }
        self.set += 1;
        for (v, o) in self.v.iter_mut().zip(o.v.iter()) {
            v.v = v.v.max(o.v);
        }
    }

    pub fn min(&mut self, o: &PeaAssVar) {
        if self.set == 0 {
            for (v, o) in self.v.iter_mut().zip(o.v.iter()) {
                v.v = o.v;
            }
        }
        self.set += 1;
        for (v, o) in self.v.iter_mut().zip(o.v.iter()) {
            v.v = v.v.max(o.v);
        }
    }
    pub fn weigh(&mut self, o: &PeaAssVar) {
        for (v, o) in self.v.iter_mut().zip(o.v.iter()) {
            v.v *= o.v;
        }
    }
    pub fn sum(&mut self) {
        self.res = self.v.iter().map(|v| v.v).sum();
    }
}
use crate::p01::AojObj;
impl PeaTrans {
    pub fn from_cmt(&mut self, cmt: &sglib04::geo1::CnlData) {
        self.tr_tag = cmt.tr_tag.clone();
        self.tr_fid = cmt.tr_fid.clone();
        self.tr_lt = cmt.tr_lt;
        self.tr_ln = cmt.tr_ln;
        self.tr_cd = cmt.tr_cd;
        self.tr_aoj = cmt.tr_aoj.clone();
        self.tr_pea = cmt.tr_pea.clone();
        self.tr_kva = cmt.tr_kva;
        self.tr_own = cmt.tr_own.clone();
        self.tr_loc = cmt.tr_loc.clone();
        self.tr_n1d = cmt.tr_n1d;
    }
}

impl PeaMeter {
    pub fn from_cmt(&mut self, cmt: &sglib04::geo1::CnlData) {
        self.mt_ins = cmt.mt_ins.clone();
        self.mt_pea = cmt.mt_pea.clone();
        self.mt_tag = cmt.mt_tag.clone();
        self.mt_phs = cmt.mt_phs.clone();
        self.mt_x = cmt.mt_x;
        self.mt_y = cmt.mt_y;
        self.mt_lt = cmt.mt_lt;
        self.mt_ln = cmt.mt_ln;
        self.mt_aoj = cmt.mt_aoj.clone();
        self.tr_tag = cmt.tr_tag.clone();
        self.tr_fid = cmt.tr_fid.clone();
        self.tr_lt = cmt.tr_lt;
        self.tr_ln = cmt.tr_ln;
        self.tr_cd = cmt.tr_cd;
        self.tr_aoj = cmt.tr_aoj.clone();
        self.tr_pea = cmt.tr_pea.clone();
        self.tr_kva = cmt.tr_kva;
        self.tr_own = cmt.tr_own.clone();
        self.tr_loc = cmt.tr_loc.clone();
        self.tr_n1d = cmt.tr_n1d;
        self.mt_n1d = cmt.mt_n1d;
        self.ar = cmt.ar.clone();
        self.ly = cmt.ly.clone();
        self.ix = cmt.ix;
    }
    pub fn from_bil(&mut self, bil: &sglib04::geo1::MeterBill) {
        self.trsg = bil.trsg.clone();
        self.pea = bil.pea.clone();
        self.ca = bil.ca.clone();
        self.inst = bil.inst.clone();
        self.rate = bil.rate.clone();
        self.volt = bil.volt.clone();
        self.mru = bil.mru.clone();
        self.mat = bil.mat.clone();
        self.main = bil.main.clone();
        self.kwh15 = bil.kwh15;
        self.kwh18 = bil.kwh18;
        self.amt19 = bil.amt19;
        self.ar = bil.ar.clone();
        self.idx = bil.idx;
        self.meth = bil.meth;
    }
}

use crate::p01::get_tr_volta;
use crate::p01::mon_kwh_2_kw;
use crate::p01::trf_kva_2_kw;
use sglib04::geo4::PowerProdType;
use sglib04::geo4::GPPS;

///
/// g0 = ProcEngine::prep1();
///
/// evpv: p13_ev_distr(&EV_PRV_ADJ_1)
///
/// sbif: sub_inf()
///
/// lp23: p03_load_lp("2023")
///
/// lp24: p03_load_lp("2024")
///
/// for (sb, sf) in &g0.sbif {
///
///   for id in &aids {
///
/// dnm = "/mnt/e/CHMBACK/pea-data/c01_pea";
///
/// bin: `Vec<u8>` = bincode::encode_to_vec(&pea, bincode::config::standard())
///
/// write(format!("{dnm}/000_pea.bin"), bin)
///
/// bin: `Vec<u8>` = bincode::encode_to_vec(&sb, bincode::config::standard())
///
/// std::fs::write(format!("{dnm}/{}.bin", sb.sbid), bin)
///
pub fn c01_chk_01() -> Result<(), Box<dyn Error>> {
    let dnm = "/mnt/e/CHMBACK/pea-data/c01_pea";
    let g0 = ProcEngine::prep5();
    let pea = c01_chk_01_01(dnm, &g0)?;
    c01_chk_01_02(&pea, dnm, &g0)?;
    Ok(())
}

pub fn c01_chk_01_01(dnm: &str, g0: &ProcEngine) -> Result<Pea, Box<dyn Error>> {
    std::fs::create_dir_all(dnm)?;
    let mut pea = Pea::default();
    for (sb, sf) in &g0.sbif {
        let ar = sf.arid.to_string();
        let ar_e = pea.aream.entry(ar).or_insert_with(|| PeaArea {
            arid: sf.arid.to_string(),
            ..Default::default()
        });
        let pv_e = ar_e
            .provm
            .entry(sf.prov.to_string())
            .or_insert_with(|| PeaProv {
                pvnm: sf.prov.to_string(),
                ..Default::default()
            });

        let Some(ev) = g0.evpv.get(&sf.prov) else {
            continue;
        };
        let Some(gpp) = GPPS.get(&sf.prov) else {
            continue;
        };
        pv_e.evpc = ev.ev_pc;
        pv_e.gppv = *gpp as f32;

        let _sb_e = pv_e.subm.entry(sb.to_string()).or_insert_with(|| PeaSub {
            sbid: sb.to_string(),
            name: sf.name.clone(),
            enam: sf.enam.clone(),
            area: sf.area.clone(),
            arid: sf.arid.clone(),
            volt: sf.volt.clone(),
            cate: sf.cate.clone(),
            egat: sf.egat.clone(),
            state: sf.state.clone(),
            conf: sf.conf.clone(),
            trax: sf.trax.clone(),
            mvax: sf.mvax.clone(),
            feed: sf.feed.clone(),
            feno: sf.feno,
            feeders: sf.feeders.clone(),
            trxn: sf.trxn,
            mvxn: sf.mvxn,
            prov: sf.prov.clone(),
            ..Default::default()
        });
    }

    let bin: Vec<u8> = bincode::encode_to_vec(&pea, bincode::config::standard()).unwrap();
    std::fs::write(format!("{dnm}/000_pea.bin"), bin).unwrap();
    println!("write to 000_pea.bin");

    Ok(pea)
}

pub fn c01_chk_01_02(pea: &Pea, dnm: &str, g0: &ProcEngine) -> Result<(), Box<dyn Error>> {
    let smrt = Regex::new(r"[12].*").unwrap();
    let mut aids: Vec<_> = pea.aream.keys().collect();
    aids.sort();

    for id in &aids {
        println!("ar:{id}");
        let id = id.to_string();
        let eg = ProcEngine::prep2(&id);
        let Some(ar) = pea.aream.get(&id) else {
            continue;
        };
        let mut pids: Vec<_> = ar.provm.keys().collect();
        pids.sort();
        for id in &pids {
            let pid = id.to_string();
            let Some(pr) = ar.provm.get(&pid) else {
                continue;
            };
            //////////////////////////////////////////////
            //////////////////////////////////////////////
            //////////////////////////////////////////////
            // province
            // Car registration
            let Some(ev) = g0.evpv.get(&pid) else {
                continue;
            };
            let Some(gpp) = GPPS.get(&pid) else {
                continue;
            };
            println!("  p:{id}");
            println!("    ev rt: {}", ev.ev_pc);
            println!("    gpp : {}", gpp);

            let mut sids: Vec<_> = pr.subm.keys().collect();
            sids.sort();
            for id in &sids {
                let sbid = id.to_string();
                let sid = id.to_string();
                let Some(sb) = pr.subm.get(&sid) else {
                    println!("          ====================== NO1 {sid}");
                    continue;
                };
                let is = eg
                    .subs
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| r.sbid == sbid)
                    .map(|(i, _)| i)
                    .collect::<Vec<_>>();
                //println!("    s:{id} is:{is:?}");
                if is.is_empty() {
                    println!("          ====================== NO2 {sid}");
                    continue;
                }
                let si = is[0];
                let sub = &eg.subs[si];

                ///////////////////////////////////////////////
                //////////////////////////////////////////////
                //  Substation Info
                let mut sb = sb.clone();
                sb.sbtp = sub.conf.to_string();
                sb.n1d_s = sub.n1d_s;
                sb.n1d_f = sub.n1d_f;
                // substation info
                if let Some(slp) = g0.lp23.get(&sbid) {
                    sb.lp_rep_23 = slp.clone();
                }
                if let Some(slp) = g0.lp24.get(&sbid) {
                    sb.lp_rep_24 = slp.clone();
                }

                ////// collect VSPP under substation
                let mut vspps = vec![];
                let vsp = &eg.vssb[si];
                if !vsp.is_empty() {
                    for pi in vsp {
                        vspps.push(*pi);
                    }
                }
                sb.vspps = vspps;

                let mut spps = vec![];
                let spp = &eg.spsb[si];
                if !spp.is_empty() {
                    for pi in spp {
                        spps.push(*pi);
                    }
                }
                sb.spps = spps;

                let mut repls = vec![];
                let repl = &eg.resb[si];
                if !repl.is_empty() {
                    for pi in repl {
                        repls.push(*pi);
                    }
                }
                sb.repls = repls;

                println!("     ID:{}", sub.sbid);
                let mut fds = sub.feed.keys().collect::<Vec<_>>();
                fds.sort();
                let mut aoj_tr = HashMap::<usize, usize>::new();
                for fid in fds {
                    let fid = fid.to_string();
                    //println!("      {fid}");
                    let Some(tis) = sub.feed.get(&fid) else {
                        continue;
                    };
                    let fd = sb.feedm.entry(fid.to_string()).or_insert_with(|| PeaFeed {
                        fdid: fid.to_string(),
                        ..Default::default()
                    });
                    for ti in tis {
                        let tr = &eg.ctrs[*ti];
                        let t1d = tr.n1d;
                        let trs = fd.tranm.entry(t1d).or_insert_with(|| PeaTrans {
                            trid: tr.trid.to_string(),
                            pea: tr.pea.to_string(),
                            n1d: tr.n1d,
                            n1d_f: tr.n1d_f,
                            ix: tr.ix,
                            lix: tr.lix,
                            mts: tr.mts.clone(),
                            ..Default::default()
                        });

                        ////////////// AojObj
                        trs.aojs = eg.aotr[*ti].clone();
                        for ai in &trs.aojs {
                            let ai = *ai;
                            if let Some(cn) = aoj_tr.get_mut(&ai) {
                                *cn += 1;
                            } else {
                                aoj_tr.insert(ai, 1);
                            }
                        }
                        trs.amps = eg.amtr[*ti].clone();
                        trs.muns = eg.mutr[*ti].clone();
                        trs.zons = eg.zntr[*ti].clone();
                        trs.sols = if eg.sotr.len() > *ti {
                            eg.sotr[*ti].clone()
                        } else {
                            Vec::<_>::new()
                        };

                        let tcm = &eg.cmts[tr.ix];
                        let Some(ow) = &tcm.tr_own else {
                            continue;
                        };
                        trs.own = ow.clone();
                        trs.vols = eg.votr[*ti].clone();
                        let (vopw, vose) = get_tr_volta(*ti, &eg);
                        trs.vopw = vopw;
                        trs.vose = vose;
                        trs.from_cmt(tcm);
                        //println!("        trs: {}", trs.n1d);

                        let Some(kv) = &tcm.tr_kva else {
                            continue;
                        };
                        if *kv == 0.0 {
                            continue;
                        }
                        trs.kw = trf_kva_2_kw(*kv);

                        for mi in &tr.mts {
                            let mt = &eg.cmts[*mi];
                            let mb = &eg.m2bs[*mi];
                            if mb.is_empty() {
                                continue;
                            }
                            let bl = &eg.bils[mb[0]];
                            let tp = if smrt.captures(bl.rate.as_str()).is_some() {
                                MeterAccType::Small
                            } else {
                                MeterAccType::Large
                            };
                            let mut met = PeaMeter::default();
                            met.from_cmt(mt);
                            met.from_bil(bl);
                            met.met_type = tp;
                            trs.mets.push(met);
                        }
                    } //end of trans
                } // end of feeder
                let mut aojs: Vec<(usize, usize)> =
                    aoj_tr.into_iter().map(|(k, v)| (v, k)).collect();
                aojs.sort_by(|a, b| b.0.cmp(&a.0));
                let mut aojv = Vec::<AojObj>::new();
                for (v, ai) in aojs {
                    let ao = &eg.aojs[ai];
                    let code = ao.code.clone().unwrap_or("".to_string());
                    let sht_name = ao.sht_name.clone().unwrap_or("".to_string());
                    let office = ao.office.clone().unwrap_or("".to_string());
                    let pea = ao.pea.clone().unwrap_or("".to_string());
                    let aoj_sz = ao.aoj_sz.clone().unwrap_or("".to_string());
                    let reg = ao.reg.clone().unwrap_or("".to_string());
                    let name = ao.name.clone().unwrap_or("".to_string());
                    let level = ao.level.unwrap_or(0f32);
                    let trcn = v;
                    //aojv.push((eg.aojs[ai].clone(), v));
                    let aoj = AojObj {
                        code: code.to_string(),
                        sht_name: sht_name.to_string(),
                        office: office.to_string(),
                        pea: pea.to_string(),
                        aoj_sz: aoj_sz.to_string(),
                        reg: reg.to_string(),
                        name: name.to_string(),
                        level,
                        trcn,
                    };
                    aojv.push(aoj);
                } // end aoj loop
                sb.aojv = aojv;

                let bin: Vec<u8> =
                    bincode::encode_to_vec(&sb, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{}.bin", sb.sbid), bin).unwrap();
                println!("       ===== write to {}.bin", sb.sbid);
            } // end sub loop
        } // end province loop
    } // end area loop

    Ok(())
}

impl PeaAssVar {
    pub fn from(n1d: u64) -> Self {
        let mut v = Vec::<AssVar>::new();
        for vt in VarType::iter() {
            let st = match vt {
                VarType::MaxPosPowSubVs01 => SumType::Max,
                VarType::MaxNegPowSubVs02 => SumType::Max,
                VarType::MaxPosPowFeederVf01 => SumType::Max,
                VarType::MaxNegPowFeederVf02 => SumType::Max,
                VarType::MaxPosDiffFeederVf03 => SumType::Max,
                VarType::MaxNegDiffFeederVf04 => SumType::Max,
                VarType::UnbalPowRate => SumType::Max,
                _ => SumType::Sum,
            };
            v.push(AssVar::new(vt, st));
        }
        PeaAssVar {
            n1d,
            v,
            ..Default::default()
        }
    }
}

/// read 000_pea.bin
/// read SSS.bin
/// write
pub fn c01_chk_02() -> Result<(), Box<dyn Error>> {
    let dnm = "/mnt/e/CHMBACK/pea-data/c01_pea";
    let buf = std::fs::read(format!("{dnm}/000_pea.bin")).unwrap();
    let (pea, _): (Pea, usize) =
        bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
    println!("pea ar:{}", pea.aream.len());
    let mut aids: Vec<_> = pea.aream.keys().collect();
    aids.sort();
    //let mut tras_mx1 = PeaAssVar::default();
    //let mut tras_mx2 = PeaAssVar::default();
    let mut tras_mx1 = PeaAssVar::from(0u64);
    let mut tras_mx2 = PeaAssVar::from(0u64);
    let mut tras_sm2 = PeaAssVar::from(0u64);
    chk_02_1(&aids, &pea, dnm, &mut tras_mx1)?;
    chk_02_2(&aids, &pea, dnm, &tras_mx1, &mut tras_mx2, &mut tras_sm2)?;
    chk_02_3(&aids, &pea, dnm, &tras_mx2, &tras_sm2)?;
    Ok(())
}

use iterstats::Iterstats;

pub fn c04_chk_lp_01() -> Result<(), Box<dyn Error>> {
    let dnm = "/mnt/e/CHMBACK/pea-data/c01_pea";
    let buf = std::fs::read(format!("{dnm}/000_pea.bin")).unwrap();
    let (pea, _): (Pea, usize) =
        bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
    println!("pea ar:{}", pea.aream.len());
    let mut aids: Vec<_> = pea.aream.keys().collect();
    aids.sort();
    let e0 = ProcEngine::prep5();
    let keys: Vec<_> = e0.lp24.keys().collect();
    let re = Regex::new(r"([A-Z]{3})-([0-9]{2})[VW].*").unwrap();
    let mut fd2fd = HashMap::<String, String>::new();
    for k in keys {
        for cap in re.captures_iter(k) {
            let a = &cap[1].to_string();
            let b = &cap[2].to_string();
            let fd = format!("{a}{b}");
            if let Some(o) = fd2fd.get(&fd) {
                println!("DUP {o} => fd:{fd} k:{k}");
            } else {
                fd2fd.insert(fd, k.to_string());
            }
        }
    }
    for id in aids {
        let aid = id.to_string();
        let Some(ar) = pea.aream.get(&aid) else {
            continue;
        };
        println!("ar:{}", ar.arid);
        let mut pids: Vec<_> = ar.provm.keys().collect();
        pids.sort();
        // province loop1
        for pid in pids {
            let Some(prov) = ar.provm.get(pid) else {
                continue;
            };
            println!("  pv:{pid}");
            let mut sids: Vec<_> = prov.subm.keys().collect();
            sids.sort();
            for sid in sids {
                let Some(_sb) = prov.subm.get(sid) else {
                    continue;
                };
                /*
                let Ok(buf) = std::fs::read(format!("{dnm}/{sid}.bin")) else {
                    continue;
                };
                */
                //let (sub, _): (PeaSub, usize) =
                //    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();

                if let Some(lp) = e0.lp24.get(sid) {
                    if let Some(v) = lp.pos_rep.val {
                        let g = v.into_iter().filter_map(|t| t).collect::<Vec<_>>();
                        let m = g.iter().mean();
                        let s = g.iter().stddev();
                        let x = g
                            .iter()
                            .filter(|v| {
                                let z = (*v - m) / s;
                                z > 2.5f32
                            })
                            //.map(|&v| v)
                            .collect::<Vec<_>>();
                        if let (Some(a), Some(n)) = (
                            x.clone()
                                .into_iter()
                                .max_by(|x, y| x.partial_cmp(y).unwrap()),
                            x.clone()
                                .into_iter()
                                .min_by(|x, y| x.partial_cmp(y).unwrap()),
                        ) {
                            let c = x.len();
                            println!("      sub: 24p {sid} z:{a}-{n} c:{c} m:{m} s:{s}");
                        }
                    }
                }
                /*
                if let Some(lp) = e0.lp23.get(sid) {
                    for v in lp.pos_rep.val.into_iter().flatten() {
                        v23p = v23p.max(v.unwrap_or(0f32));
                    }
                };
                let (mut v24n, mut v23n) = (0f32, 0f32);
                if let Some(lp) = e0.lp24.get(sid) {
                    for v in lp.neg_rep.val.into_iter().flatten() {
                        v24n = v24n.max(v.unwrap_or(0f32));
                    }
                }
                if let Some(lp) = e0.lp23.get(sid) {
                    for v in lp.neg_rep.val.into_iter().flatten() {
                        v23n = v23n.max(v.unwrap_or(0f32));
                    }
                };
                */
            }
        }
    }
    Ok(())
}

use crate::p01::p01_chk;
use crate::p08::p08_class_val;
use crate::p08::ProfType;

pub fn chk_02_1(
    aids: &Vec<&String>,
    pea: &Pea,
    dnm: &str,
    tras_mx1: &mut PeaAssVar,
) -> Result<(), Box<dyn Error>> {
    let e0 = ProcEngine::prep5();
    let keys: Vec<_> = e0.lp24.keys().collect();
    let re = Regex::new(r"([A-Z]{3})-([0-9]{2})[VW].*").unwrap();
    let mut fd2fd = HashMap::<String, String>::new();
    for k in keys {
        for cap in re.captures_iter(k) {
            let a = &cap[1].to_string();
            let b = &cap[2].to_string();
            let fd = format!("{a}{b}");
            if let Some(o) = fd2fd.get(&fd) {
                println!("DUP {o} => fd:{fd} k:{k}");
            } else {
                fd2fd.insert(fd, k.to_string());
            }
        }
    }

    // area loop 1
    for id in aids {
        let aid = id.to_string();
        let Some(ar) = pea.aream.get(&aid) else {
            continue;
        };
        println!("ar:{}", ar.arid);
        let eg = ProcEngine::prep3(id);

        let mut am_dn = HashMap::<String, (f32, f32)>::new();
        let mut mu_dn = HashMap::<String, (f32, f32)>::new();
        for dn in &eg.amps {
            let key = dn.key.to_string();
            if let Some((_po, aa)) = am_dn.get_mut(&key) {
                *aa += dn.area;
            } else {
                am_dn.insert(key, (dn.popu, dn.area));
            }
        }
        for dn in &eg.muni {
            let key = dn.key.to_string();
            if let Some((_po, aa)) = mu_dn.get_mut(&key) {
                *aa += dn.area;
            } else {
                mu_dn.insert(key, (dn.popu, dn.area));
            }
        }

        let mut pids: Vec<_> = ar.provm.keys().collect();
        pids.sort();
        // province loop1
        for pid in pids {
            let Some(prov) = ar.provm.get(pid) else {
                continue;
            };
            let vp01 = prov.evpc;
            let vp02 = prov.gppv;
            let evlk = *EV_LIKELY.get(pid).unwrap_or(&0f32);
            let selk = *SELE_LIKELY.get(pid).unwrap_or(&0f32);

            println!("  pv:{pid}");
            let mut sids: Vec<_> = prov.subm.keys().collect();
            sids.sort();
            for sid in sids {
                let Some(_sb) = prov.subm.get(sid) else {
                    continue;
                };
                let Ok(buf) = std::fs::read(format!("{dnm}/{sid}.bin")) else {
                    continue;
                };
                let (sub, _): (PeaSub, usize) =
                    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
                println!("      feed: {}", sub.feeders.len());

                // Substation
                //=====================================================
                //=====================================================
                let mut vs01 = 0f32;
                if let Some(lp) = e0.lp24.get(sid) {
                    for v in lp.pos_rep.val.into_iter().flatten() {
                        vs01 = vs01.max(v.unwrap_or(0f32));
                    }
                } else if let Some(lp) = e0.lp23.get(sid) {
                    for v in lp.pos_rep.val.into_iter().flatten() {
                        vs01 = vs01.max(v.unwrap_or(0f32));
                    }
                };
                let mut vs02 = 0f32;
                if let Some(lp) = e0.lp24.get(sid) {
                    for v in lp.neg_rep.val.into_iter().flatten() {
                        vs02 = vs02.max(v.unwrap_or(0f32));
                    }
                } else if let Some(lp) = e0.lp23.get(sid) {
                    for v in lp.neg_rep.val.into_iter().flatten() {
                        vs02 = vs02.max(v.unwrap_or(0f32));
                    }
                };
                // LP24 - Phase
                ////////////////////////////
                let mut solar = 0f32;
                if let Some(lp) = e0.lp24.get(sid) {
                    if let Some(lp) = &lp.pos_rep.val {
                        if let Ok(lpf) = p08_class_val(lp) {
                            if lpf.lp_type == ProfType::SolarPower
                                || lpf.lp_type == ProfType::SolarNight
                            {
                                solar = lpf.sol_en.unwrap_or(0f32);
                            }
                        }
                        //
                        //
                    }
                }
                // LP23 - Phase
                ////////////////////////////

                //  VSPP, SPP, RE plan
                let mut vs03 = 0f32;
                for pi in &sub.vspps {
                    vs03 += eg.vsps[*pi].kw.unwrap_or(0f32);
                }
                let mut vs04 = 0f32;
                for pi in &sub.spps {
                    vs04 += eg.spps[*pi].mw.unwrap_or(0f32);
                }
                let mut vs05 = 0f32;
                for pi in &sub.repls {
                    if let PowerProdType::SPP = eg.repl[*pi].pptp {
                        vs05 += eg.repl[*pi].pwmw.unwrap_or(0f32);
                    }
                }
                let mut vs06 = 0f32;
                for pi in &sub.repls {
                    if let PowerProdType::VSPP = eg.repl[*pi].pptp {
                        vs06 += eg.repl[*pi].pwmw.unwrap_or(0f32);
                    }
                }
                let vs07 = sub.mvxn as f32;

                println!("    sb:{sid} feed: {}", sub.feeders.len());
                let mut fids: Vec<_> = sub.feedm.keys().collect();
                fids.sort();
                let mut s_tr_ass = Vec::<PeaAssVar>::new();
                for fid in fids {
                    let Some(fd) = sub.feedm.get(fid) else {
                        continue;
                    };
                    ////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////
                    // Feeder
                    let (mut af01, mut af03, mut af02, mut af04) = (None, None, None, None);
                    let k1 = format!("{fid}");
                    let key = fd2fd.get(&k1).unwrap_or(&"-".to_string()).to_string();
                    if let Some(lp) = e0.lp24.get(&key) {
                        if let Some(vv) = lp.pos_rep.val {
                            for v in vv.iter().flatten() {
                                if let Some(v0) = af01 {
                                    af01 = Some(v.max(v0))
                                } else {
                                    af01 = Some(*v);
                                }
                                if let Some(v0) = af03 {
                                    af03 = Some(v.min(v0))
                                } else {
                                    af03 = Some(*v);
                                }
                            }
                        }
                        if let Some(vv) = lp.neg_rep.val {
                            for v in vv.iter().flatten() {
                                if let Some(v0) = af02 {
                                    af02 = Some(v.max(v0))
                                } else {
                                    af02 = Some(*v);
                                }
                                if let Some(v0) = af04 {
                                    af04 = Some(v.min(v0))
                                } else {
                                    af04 = Some(*v);
                                }
                            }
                        }
                    } else if let Some(lp) = e0.lp23.get(&key) {
                        if let Some(vv) = lp.pos_rep.val {
                            for v in vv.iter().flatten() {
                                if let Some(v0) = af01 {
                                    af01 = Some(v.max(v0))
                                } else {
                                    af01 = Some(*v);
                                }
                                if let Some(v0) = af03 {
                                    af03 = Some(v.min(v0))
                                } else {
                                    af03 = Some(*v);
                                }
                            }
                        }
                        if let Some(vv) = lp.neg_rep.val {
                            for v in vv.iter().flatten() {
                                if let Some(v0) = af02 {
                                    af02 = Some(v.max(v0))
                                } else {
                                    af02 = Some(*v);
                                }
                                if let Some(v0) = af04 {
                                    af04 = Some(v.min(v0))
                                } else {
                                    af04 = Some(*v);
                                }
                            }
                        }
                    };
                    let vf01 = af01.unwrap_or(0f32);
                    let mut vf03 = af03.unwrap_or(0f32);
                    vf03 = vf01 - vf03;
                    let vf02 = af02.unwrap_or(0f32);
                    let mut vf04 = af04.unwrap_or(0f32);
                    vf04 = vf02 - vf04;

                    let mut tids: Vec<_> = fd.tranm.keys().collect();
                    tids.sort();
                    for tid in tids {
                        let Some(trn) = fd.tranm.get(tid) else {
                            continue;
                        };
                        let aojs = trn.aojs.clone();
                        let vt05 = trn.tr_kva.unwrap_or(10f32);
                        let vt05 = trf_kva_2_kw(vt05);
                        let mut vt06 = 1f32;
                        for zi in &trn.zons {
                            match eg.zons[*zi].zncd.clone().expect("-").as_str() {
                                "21" | "22" | "24" => {
                                    vt06 = vt06.max(5f32);
                                }
                                "11" | "12" | "13" | "14" => {
                                    vt06 = vt06.max(4f32);
                                }
                                "23" | "25" | "31" => {
                                    vt06 = vt06.max(3f32);
                                }
                                "41" | "42" => {
                                    vt06 = vt06.max(2f32);
                                }
                                _ => {}
                            }
                        }

                        let aoj = if aojs.is_empty() {
                            "-".to_string()
                        } else {
                            let ai = aojs[0];
                            eg.aojs[ai]
                                .sht_name
                                .clone()
                                .unwrap_or("-".to_string())
                                .to_string()
                        };
                        let mut vt01 = 0f32;
                        let mut vt02 = 0f32;
                        let mut vt10 = 0f32;

                        let (mut se_a, mut se_b, mut se_c) = (0.0, 0.0, 0.0);
                        let (mut sl_a, mut sl_b, mut sl_c, mut sl_3) = (0.0, 0.0, 0.0, 0.0);
                        for met in &trn.mets {
                            ///////////////////////////////////////////////////
                            ///////////////////////////////////////////////////
                            // Meter
                            if let MeterAccType::Small = met.met_type {
                                if met.main.is_empty() && met.kwh18 > 200f32 {
                                    vt01 += 1f32;
                                    vt02 += met.kwh15;
                                }
                            } else if let MeterAccType::Large = met.met_type {
                                vt10 += met.kwh15;
                                print!("_{}", met.kwh15);
                            }
                            if trn.own == "P" {
                                match met.mt_phs.clone().unwrap_or(String::new()).as_str() {
                                    "A" => se_a += met.kwh15,
                                    "B" => se_b += met.kwh15,
                                    "C" => se_c += met.kwh15,
                                    _ => {}
                                }
                                match met.mt_phs.clone().unwrap_or(String::new()).as_str() {
                                    "A" => sl_a += met.kwh15,
                                    "B" => sl_b += met.kwh15,
                                    "C" => sl_c += met.kwh15,
                                    _ => sl_3 += met.kwh15,
                                }
                            }
                        }
                        let vt11 = trn.mets.len() as f32;
                        let vt12 = 1f32;
                        sl_3 = mon_kwh_2_kw(sl_3);
                        sl_a = mon_kwh_2_kw(sl_a);
                        sl_b = mon_kwh_2_kw(sl_b);
                        sl_c = mon_kwh_2_kw(sl_c);
                        let v_phs_a = sl_3 + sl_a;
                        let v_phs_b = sl_3 + sl_b;
                        let v_phs_c = sl_3 + sl_c;
                        let v_all_p = sl_3 + sl_a + sl_b + sl_c;
                        let v_ph_av = (v_phs_a + v_phs_b + v_phs_c) / 3f32;
                        let v_ph_mx = v_phs_a.max(v_phs_b.max(v_phs_c));
                        let v_ph_rt = v_ph_mx / z2o(v_ph_av);
                        let v_al_kw = v_phs_a + v_phs_b + v_phs_c;
                        const TRF_LOSS_RATIO: f32 = 0.03;
                        const TRF_UNBAL_K: f32 = 1.0f32;
                        const TRF_UNBAL_CNT_RATE: f32 = 0.8f32;
                        let v_loss = v_al_kw * TRF_LOSS_RATIO;
                        let v_unba = v_loss * TRF_UNBAL_K * v_ph_rt * v_ph_rt;
                        let v_unb_sat = v_ph_mx / z2o(vt05);
                        let v_unb_cnt = if v_unb_sat >= TRF_UNBAL_CNT_RATE {
                            1f32
                        } else {
                            0f32
                        };
                        let v_max_sat = v_all_p / z2o(vt05);
                        let v_max_cnt = if v_unb_cnt == 0f32 && v_max_sat >= TRF_UNBAL_CNT_RATE {
                            1f32
                        } else {
                            0f32
                        };

                        let mut vt08 = 0f32;
                        let se_p = se_a + se_b + se_c;
                        if se_a < se_p && se_b < se_p && se_c < se_p {
                            let ab = (se_a - se_b).abs();
                            let bc = (se_b - se_c).abs();
                            let ca = (se_c - se_a).abs();
                            vt08 = (ab + bc + ca) * 0.5;
                        }
                        let vt08 = mon_kwh_2_kw(vt08);
                        //let vt09 = trf_kva_2_kw(vt02);
                        let vt09 = mon_kwh_2_kw(vt02);

                        let mut vt03 = 0f32;
                        for vi in &trn.vols {
                            for (pw, no) in &eg.vols[*vi].chgr {
                                vt03 += (*pw * *no) as f32;
                            }
                        }
                        let mut vt04 = 0f32;
                        for vi in &trn.vols {
                            for (_yr, am) in &eg.vols[*vi].sell {
                                vt04 += *am;
                            }
                        }
                        let mut vt07 = 1f32;
                        for ai in &trn.amps {
                            let am = &eg.amps[*ai].key;
                            if let Some((p, a)) = am_dn.get(am) {
                                let a = a / 1_000f32;
                                let pd = p / a * 0.6f32;
                                let v = match pd {
                                    0f32..30f32 => 1f32,
                                    30f32..60f32 => 2f32,
                                    60f32..150f32 => 3f32,
                                    150f32..500f32 => 4f32,
                                    _ => 5f32,
                                };
                                vt07 = vt07.max(v);
                            }
                        }
                        for ai in &trn.muns {
                            let mu = &eg.muni[*ai].key;
                            if let Some((p, a)) = mu_dn.get(mu) {
                                let a = a / 1_000f32;
                                let pd = p / a * 2.5f32;
                                let v = match pd {
                                    0f32..15f32 => 6f32,
                                    15f32..30f32 => 7f32,
                                    30f32..70f32 => 8f32,
                                    70f32..200f32 => 9f32,
                                    _ => 10f32,
                                };
                                vt07 = vt07.max(v);
                            }
                        }

                        let mut tr_as = PeaAssVar::from(trn.n1d);
                        tr_as.arid = aid.to_string();
                        tr_as.pvid = pid.to_string();
                        tr_as.sbid = sid.to_string();
                        tr_as.fdid = fid.to_string();
                        tr_as.own = trn.own.to_string();
                        tr_as.peano = trn.tr_pea.clone().unwrap_or("".to_string()).to_string();
                        tr_as.aoj = aoj;
                        tr_as.v[VarType::None as usize].v = 0f32;
                        tr_as.v[VarType::NewCarRegVp01 as usize].v = vp01;
                        tr_as.v[VarType::GppVp02 as usize].v = vp02;
                        tr_as.v[VarType::MaxPosPowSubVs01 as usize].v = vs01;
                        tr_as.v[VarType::MaxNegPowSubVs02 as usize].v = vs02;
                        tr_as.v[VarType::VsppMvVs03 as usize].v = vs03;
                        tr_as.v[VarType::SppHvVs04 as usize].v = vs04;
                        tr_as.v[VarType::BigLotMvVs05 as usize].v = vs05;
                        tr_as.v[VarType::BigLotHvVs06 as usize].v = vs06;
                        tr_as.v[VarType::SubPowCapVs07 as usize].v = vs07;
                        tr_as.v[VarType::MaxPosPowFeederVf01 as usize].v = vf01;
                        tr_as.v[VarType::MaxNegPowFeederVf02 as usize].v = vf02;
                        tr_as.v[VarType::MaxPosDiffFeederVf03 as usize].v = vf03;
                        tr_as.v[VarType::MaxNegDiffFeederVf04 as usize].v = vf04;
                        tr_as.v[VarType::NoMeterTransVt01 as usize].v = vt01;
                        tr_as.v[VarType::SmallSellTrVt02 as usize].v = vt02;
                        tr_as.v[VarType::ChgStnCapTrVt03 as usize].v = vt03;
                        tr_as.v[VarType::ChgStnSellTrVt04 as usize].v = vt04;
                        tr_as.v[VarType::PwCapTriVt05 as usize].v = vt05;
                        tr_as.v[VarType::ZoneTrVt06 as usize].v = vt06;
                        tr_as.v[VarType::PopTrVt07 as usize].v = vt07;
                        tr_as.v[VarType::UnbalPowTrVt08 as usize].v = vt08;
                        tr_as.v[VarType::PkPowTrVt09 as usize].v = vt09;
                        tr_as.v[VarType::LargeSellTrVt10 as usize].v = vt10;
                        tr_as.v[VarType::AllNoMeterTrVt11 as usize].v = vt11;
                        tr_as.v[VarType::NoTrVt12 as usize].v = vt12;
                        tr_as.v[VarType::PkSelPowPhsAKw as usize].v = v_phs_a;
                        tr_as.v[VarType::PkSelPowPhsBKw as usize].v = v_phs_b;
                        tr_as.v[VarType::PkSelPowPhsCKw as usize].v = v_phs_c;
                        tr_as.v[VarType::PkSelPowPhsAvg as usize].v = v_ph_av;
                        tr_as.v[VarType::PkSelPowPhsMax as usize].v = v_ph_mx;
                        tr_as.v[VarType::UnbalPowRate as usize].v = v_ph_rt;
                        tr_as.v[VarType::TransLossKw as usize].v = v_loss;
                        tr_as.v[VarType::UnbalPowLossKw as usize].v = v_unba;
                        tr_as.v[VarType::CntTrUnbalLoss as usize].v = v_unb_cnt;
                        tr_as.v[VarType::CntTrSatLoss as usize].v = v_max_cnt;
                        tr_as.v[VarType::EvCarLikely as usize].v = evlk;
                        tr_as.v[VarType::SelectLikely as usize].v = selk;
                        tr_as.v[VarType::SolarEnergy as usize].v = solar;
                        tras_mx1.max(&tr_as);
                        //s_tr_sum.add(&tr_as);
                        s_tr_ass.push(tr_as);
                    } // end trans loop
                } // end feeder loop
                write_trn_ass_01(&s_tr_ass, &format!("{dnm}/{sid}-raw0.txt"))?;
                write_ass_csv_01(&s_tr_ass, &format!("{dnm}/{sid}-raw0.csv"))?;
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&s_tr_ass, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-raw.bin"), bin).unwrap();
            } // end sub loop
        } // end provi loop
    } // end area
    Ok(())
}

pub fn chk_02_2(
    aids: &Vec<&String>,
    pea: &Pea,
    dnm: &str,
    tras_mx1: &PeaAssVar,
    tras_mx2: &mut PeaAssVar,
    tras_sm2: &mut PeaAssVar,
) -> Result<(), Box<dyn Error>> {
    //////////////////////////////////////////////
    // EV Weight
    let mut we_ev = PeaAssVar::from(0u64);
    we_ev.v[VarType::NewCarRegVp01 as usize].v = 0.20;
    we_ev.v[VarType::GppVp02 as usize].v = 0.20;
    we_ev.v[VarType::NoMeterTransVt01 as usize].v = 0.05;
    we_ev.v[VarType::SmallSellTrVt02 as usize].v = 0.30;
    we_ev.v[VarType::ChgStnCapTrVt03 as usize].v = 0.15;
    we_ev.v[VarType::ChgStnSellTrVt04 as usize].v = 0.15;
    we_ev.v[VarType::PwCapTriVt05 as usize].v = 0.00;
    we_ev.v[VarType::ZoneTrVt06 as usize].v = 0.20;
    we_ev.v[VarType::PopTrVt07 as usize].v = 0.10;
    we_ev.v[VarType::EvCarLikely as usize].v = 0.00;

    //////////////////////////////////////////////
    // Solar Weight
    let mut we_so = PeaAssVar::from(0u64);
    we_so.v[VarType::GppVp02 as usize].v = 0.20;
    /*
    we_so.v[VarType::MaxPosPowSubVs01 as usize].v = 0.02;
    we_so.v[VarType::MaxNegPowSubVs02 as usize].v = 0.02;
    we_so.v[VarType::VsppMvVs03 as usize].v = 0.02;
    we_so.v[VarType::SppHvVs04 as usize].v = 0.02;
    we_so.v[VarType::BigLotMvVs05 as usize].v = 0.02;
    we_so.v[VarType::MaxNegPowFeederVf02 as usize].v = 0.05;
    we_so.v[VarType::MaxPosDiffFeederVf03 as usize].v = 0.05;
    */
    we_so.v[VarType::NoMeterTransVt01 as usize].v = 0.05;
    we_so.v[VarType::SmallSellTrVt02 as usize].v = 0.30;
    //we_so.v[VarType::PwCapTriVt05 as usize].v = 0.20;
    we_so.v[VarType::ZoneTrVt06 as usize].v = 0.20;
    we_so.v[VarType::PopTrVt07 as usize].v = 0.10;

    for id in aids {
        let aid = id.to_string();
        let Some(ar) = pea.aream.get(&aid) else {
            continue;
        };
        println!("ar:{}", ar.arid);
        let mut pids: Vec<_> = ar.provm.keys().collect();
        pids.sort();
        for pid in pids {
            let Some(prov) = ar.provm.get(pid) else {
                continue;
            };
            println!("  pv:{pid}");
            let mut sids: Vec<_> = prov.subm.keys().collect();
            sids.sort();
            for sid in sids {
                let Some(_sb) = prov.subm.get(sid) else {
                    continue;
                };

                ////////////////////////////////////////////////
                ////////////////////////////////////////////////
                // read raw data
                let Ok(buf) = std::fs::read(format!("{dnm}/{sid}-raw.bin")) else {
                    continue;
                };
                let (mut v_tras_raw, _): (Vec<PeaAssVar>, usize) =
                    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
                println!("   {sid} - {}", v_tras_raw.len());
                // normalize data
                let mut v_tras_nor = v_tras_raw.clone();
                for tras in &mut v_tras_nor {
                    tras.nor(tras_mx1);
                }
                //// save normal bin
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras_nor, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-nor.bin"), bin).unwrap();
                //let a = write_trn_ass_01(&v_tras_nor, &format!("{dnm}/{sid}-nor.txt"))?;
                write_trn_ass_01(&v_tras_nor, &format!("{dnm}/{sid}-nor0.txt"))?;
                write_ass_csv_01(&v_tras_nor, &format!("{dnm}/{sid}-nor0.csv"))?;
                //println!("=====================================  {}", a == b);

                ////////////////////////////////////////////////
                ////////////////////////////////////////////////
                // calculate EV
                let mut v_tras_ev = v_tras_nor.clone();
                for (tras, tras0) in v_tras_ev.iter_mut().zip(v_tras_raw.iter_mut()) {
                    tras.weigh(&we_ev);
                    tras.sum();
                    //tras0.vc01 = tras.sum;
                    tras0.v[VarType::HmChgEvTrVc01 as usize].v = tras.res;
                }
                //// save ev bin
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras_ev, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-ev.bin"), bin).unwrap();
                //let a = write_trn_ass_01(&v_tras_ev, &format!("{dnm}/{sid}-ev.txt"))?;
                write_trn_ass_01(&v_tras_ev, &format!("{dnm}/{sid}-ev0.txt"))?;
                write_ass_csv_01(&v_tras_ev, &format!("{dnm}/{sid}-ev0.csv"))?;
                //println!("=====================================  {}", a == b);

                ////////////////////////////////////////////////
                ////////////////////////////////////////////////
                // calculate solar
                let mut v_tras_so = v_tras_nor.clone();
                for (tras, tras0) in v_tras_so.iter_mut().zip(v_tras_raw.iter_mut()) {
                    tras.weigh(&we_so);
                    tras.sum();
                    tras0.v[VarType::PowSolarVc07 as usize].v = tras.res;
                }
                //// save ev bin
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras_ev, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-so.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras_ev, &format!("{dnm}/{sid}-so0.txt"))?;

                ///////////////////////////////////////////////
                // summary of all data
                for tras in v_tras_raw.iter() {
                    tras_sm2.add(tras);
                }

                for tr in v_tras_raw.iter_mut() {
                    tr.v[VarType::LvPowSatTrVc02 as usize].v = tr.v[VarType::PkPowTrVt09 as usize]
                        .v
                        / z2o(tr.v[VarType::PwCapTriVt05 as usize].v);
                    tr.v[VarType::CntLvPowSatTrVc03 as usize].v =
                        if tr.v[VarType::LvPowSatTrVc02 as usize].v > 0.8f32 {
                            1f32
                        } else {
                            0f32
                        };
                    tr.v[VarType::ChgStnCapVc04 as usize].v =
                        tr.v[VarType::ChgStnCapTrVt03 as usize].v;
                    tr.v[VarType::ChgStnSellVc05 as usize].v =
                        tr.v[VarType::ChgStnSellTrVt04 as usize].v;
                    tr.v[VarType::MvPowSatTrVc06 as usize].v =
                        tr.v[VarType::MaxPosPowSubVs01 as usize].v
                            / z2o(tr.v[VarType::SubPowCapVs07 as usize].v);
                    tr.v[VarType::MvVsppVc08 as usize].v = tr.v[VarType::VsppMvVs03 as usize].v;
                    tr.v[VarType::HvSppVc09 as usize].v = tr.v[VarType::SppHvVs04 as usize].v;
                    tr.v[VarType::SmallSellVc10 as usize].v =
                        tr.v[VarType::SmallSellTrVt02 as usize].v;
                    tr.v[VarType::LargeSellVc11 as usize].v =
                        tr.v[VarType::LargeSellTrVt10 as usize].v;
                    tr.v[VarType::UnbalPowVc12 as usize].v =
                        tr.v[VarType::UnbalPowTrVt08 as usize].v;
                    let v = tr.v[VarType::UnbalPowTrVt08 as usize].v
                        / z2o(tr.v[VarType::PwCapTriVt05 as usize].v);
                    tr.v[VarType::CntUnbalPowVc13 as usize].v =
                        if v > 0.5f32 { 1f32 } else { 0f32 };

                    tras_mx2.max(tr);
                }
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras_raw, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-rw2.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras_raw, &format!("{dnm}/{sid}-rw20.txt"))?;
                write_ass_csv_01(&v_tras_raw, &format!("{dnm}/{sid}-rw20.csv"))?;
            } // end sub loop
        } // end provi loop
    } // end area
    Ok(())
}

pub fn chk_02_3(
    aids: &Vec<&String>,
    pea: &Pea,
    dnm: &str,
    tras_mx2: &PeaAssVar,
    tras_sm2: &PeaAssVar,
) -> Result<(), Box<dyn Error>> {
    let mut we_uc1 = PeaAssVar::from(0u64);
    we_uc1.v[VarType::SmallSellTrVt02 as usize].v = 0.05;
    we_uc1.v[VarType::LargeSellTrVt10 as usize].v = 0.05;
    we_uc1.v[VarType::HmChgEvTrVc01 as usize].v = 0.30;
    we_uc1.v[VarType::CntLvPowSatTrVc03 as usize].v = 0.15;
    we_uc1.v[VarType::ChgStnCapVc04 as usize].v = 0.05;
    we_uc1.v[VarType::ChgStnSellVc05 as usize].v = 0.05;
    we_uc1.v[VarType::MvPowSatTrVc06 as usize].v = 0.05;
    we_uc1.v[VarType::PowSolarVc07 as usize].v = 0.15;
    we_uc1.v[VarType::ZoneTrVt06 as usize].v = 0.05;
    we_uc1.v[VarType::PopTrVt07 as usize].v = 0.05;
    we_uc1.v[VarType::MvVsppVc08 as usize].v = 0.05;
    we_uc1.v[VarType::HvSppVc09 as usize].v = 0.05;
    we_uc1.v[VarType::UnbalPowVc12 as usize].v = 0.10;
    we_uc1.v[VarType::CntUnbalPowVc13 as usize].v = 0.05;
    we_uc1.v[VarType::SelectLikely as usize].v = 0.00;

    let mut we_uc2 = PeaAssVar::from(0u64);
    we_uc2.v[VarType::SmallSellTrVt02 as usize].v = 0.05;
    we_uc2.v[VarType::LargeSellTrVt10 as usize].v = 0.10;
    we_uc2.v[VarType::HmChgEvTrVc01 as usize].v = 0.05;
    we_uc2.v[VarType::CntLvPowSatTrVc03 as usize].v = 0.05;
    we_uc2.v[VarType::ChgStnCapVc04 as usize].v = 0.10;
    we_uc2.v[VarType::ChgStnSellVc05 as usize].v = 0.10;
    we_uc2.v[VarType::MvPowSatTrVc06 as usize].v = 0.15;
    we_uc2.v[VarType::PowSolarVc07 as usize].v = 0.05;
    we_uc2.v[VarType::MvVsppVc08 as usize].v = 0.15;
    we_uc2.v[VarType::HvSppVc09 as usize].v = 0.10;
    we_uc2.v[VarType::UnbalPowVc12 as usize].v = 0.05;
    we_uc2.v[VarType::CntUnbalPowVc13 as usize].v = 0.05;
    we_uc2.v[VarType::SelectLikely as usize].v = 0.10;

    let mut we_uc3 = PeaAssVar::from(0u64);
    we_uc3.v[VarType::SmallSellTrVt02 as usize].v = 0.10;
    we_uc3.v[VarType::LargeSellTrVt10 as usize].v = 0.05;
    we_uc3.v[VarType::HmChgEvTrVc01 as usize].v = 0.15;
    we_uc3.v[VarType::CntLvPowSatTrVc03 as usize].v = 0.15;
    we_uc3.v[VarType::ChgStnCapVc04 as usize].v = 0.05;
    we_uc3.v[VarType::ChgStnSellVc05 as usize].v = 0.05;
    we_uc3.v[VarType::MvPowSatTrVc06 as usize].v = 0.05;
    we_uc3.v[VarType::PowSolarVc07 as usize].v = 0.15;
    we_uc3.v[VarType::MvVsppVc08 as usize].v = 0.05;
    we_uc3.v[VarType::HvSppVc09 as usize].v = 0.05;
    we_uc3.v[VarType::UnbalPowVc12 as usize].v = 0.10;
    we_uc3.v[VarType::CntUnbalPowVc13 as usize].v = 0.05;
    we_uc3.v[VarType::SelectLikely as usize].v = 0.20;

    for id in aids {
        let aid = id.to_string();
        let Some(ar) = pea.aream.get(&aid) else {
            continue;
        };
        println!("ar:{}", ar.arid);
        let mut pids: Vec<_> = ar.provm.keys().collect();
        pids.sort();
        for pid in pids {
            let Some(prov) = ar.provm.get(pid) else {
                continue;
            };
            println!("  pv:{pid}");
            let mut sids: Vec<_> = prov.subm.keys().collect();
            sids.sort();
            for sid in sids {
                let Some(_sb) = prov.subm.get(sid) else {
                    continue;
                };

                ////////////////////////////////////////////////
                ////////////////////////////////////////////////
                // read raw data
                let Ok(buf) = std::fs::read(format!("{dnm}/{sid}-rw2.bin")) else {
                    continue;
                };
                let (mut v_tras_raw, _): (Vec<PeaAssVar>, usize) =
                    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
                println!("   {sid} - {}", v_tras_raw.len());
                // normalize data
                let mut v_tras_nor = v_tras_raw.clone();
                for tras in &mut v_tras_nor {
                    tras.nor(tras_mx2);
                }

                ///////////////////////////////////////////////
                // calculate ratio with the whole
                let mut v_tras_sum = v_tras_raw.clone();
                for (tras, tras0) in v_tras_sum.iter_mut().zip(v_tras_raw.iter_mut()) {
                    tras.nor(tras_sm2);
                    tras0.v[VarType::NoHmChgEvTr as usize].v =
                        tras.v[VarType::HmChgEvTrVc01 as usize].v * 210_000f32;
                    tras0.v[VarType::PowHmChgEvTr as usize].v =
                        tras0.v[VarType::NoHmChgEvTr as usize].v * 0.007f32;
                }

                //// save normal bin
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras_nor, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-no2.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras_nor, &format!("{dnm}/{sid}-no20.txt"))?;
                write_ass_csv_01(&v_tras_nor, &format!("{dnm}/{sid}-no20.csv"))?;

                //// UC1
                let mut v_uc1 = v_tras_nor.clone();
                for (tras, tras0) in v_uc1.iter_mut().zip(v_tras_raw.iter_mut()) {
                    tras.weigh(&we_uc1);
                    tras.sum();
                    //tras0.vc14 = tras.sum;
                    tras0.v[VarType::Uc1ValVc14 as usize].v = tras.res;
                }
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_uc1, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-uc1.bin"), bin).unwrap();
                write_trn_ass_01(&v_uc1, &format!("{dnm}/{sid}-uc10.txt"))?;
                write_ass_csv_01(&v_uc1, &format!("{dnm}/{sid}-uc10.csv"))?;

                //// UC2
                let mut v_uc2 = v_tras_nor.clone();
                for (tras, tras0) in v_uc2.iter_mut().zip(v_tras_raw.iter_mut()) {
                    tras.weigh(&we_uc2);
                    tras.sum();
                    //tras0.vc15 = tras.sum;
                    tras0.v[VarType::Uc2ValVc15 as usize].v = tras.res;
                }
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_uc2, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-uc2.bin"), bin).unwrap();
                write_trn_ass_01(&v_uc2, &format!("{dnm}/{sid}-uc20.txt"))?;
                write_ass_csv_01(&v_uc2, &format!("{dnm}/{sid}-uc20.csv"))?;

                //// UC3
                let mut v_uc3 = v_tras_nor.clone();
                for (tras, tras0) in v_uc3.iter_mut().zip(v_tras_raw.iter_mut()) {
                    tras.weigh(&we_uc3);
                    tras.sum();
                    //tras0.vc16 = tras.sum;
                    tras0.v[VarType::Uc3ValVc16 as usize].v = tras.res;
                }
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_uc3, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-uc3.bin"), bin).unwrap();
                write_trn_ass_01(&v_uc3, &format!("{dnm}/{sid}-uc30.txt"))?;
                write_ass_csv_01(&v_uc3, &format!("{dnm}/{sid}-uc30.csv"))?;

                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras_raw, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-rw3.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras_raw, &format!("{dnm}/{sid}-rw30.txt"))?;
                write_ass_csv_01(&v_tras_raw, &format!("{dnm}/{sid}-rw30.csv"))?;
            }
        }
    }
    Ok(())
}

fn write_trn_ass_01(tr_as: &Vec<PeaAssVar>, fnm: &str) -> Result<String, Box<dyn Error>> {
    let flds = [
        VarType::NewCarRegVp01,
        VarType::GppVp02,
        VarType::MaxPosPowSubVs01,
        VarType::MaxNegPowSubVs02,
        VarType::VsppMvVs03,
        VarType::SppHvVs04,
        VarType::BigLotMvVs05,
        VarType::BigLotHvVs06,
        VarType::SubPowCapVs07,
        VarType::MaxPosPowFeederVf01,
        VarType::MaxNegPowFeederVf02,
        VarType::MaxPosDiffFeederVf03,
        VarType::MaxNegDiffFeederVf04,
        VarType::NoMeterTransVt01,
        VarType::SmallSellTrVt02,
        VarType::ChgStnCapTrVt03,
        VarType::ChgStnSellTrVt04,
        VarType::PwCapTriVt05,
        VarType::ZoneTrVt06,
        VarType::PopTrVt07,
        VarType::UnbalPowTrVt08,
        VarType::PkPowTrVt09,
        VarType::LargeSellTrVt10,
        VarType::AllNoMeterTrVt11,
        VarType::NoTrVt12,
        VarType::HmChgEvTrVc01,
        VarType::LvPowSatTrVc02,
        VarType::CntLvPowSatTrVc03,
        VarType::ChgStnCapVc04,
        VarType::ChgStnSellVc05,
        VarType::MvPowSatTrVc06,
        VarType::PowSolarVc07,
        VarType::MvVsppVc08,
        VarType::HvSppVc09,
        VarType::SmallSellVc10,
        VarType::LargeSellVc11,
        VarType::UnbalPowVc12,
        VarType::CntUnbalPowVc13,
        VarType::Uc1ValVc14,
        VarType::Uc2ValVc15,
        VarType::Uc3ValVc16,
        VarType::Uc1RankVc17,
        VarType::Uc2RankVc18,
        VarType::Uc3RankVc19,
        VarType::NoHmChgEvTr,
        VarType::PowHmChgEvTr,
        VarType::PkSelPowPhsAKw,
        VarType::PkSelPowPhsBKw,
        VarType::PkSelPowPhsCKw,
        VarType::PkSelPowPhsAvg,
        VarType::PkSelPowPhsMax,
        VarType::UnbalPowRate,
        VarType::TransLossKw,
        VarType::UnbalPowLossKw,
    ];
    write_text_01(tr_as, &flds, fnm)
}

fn write_trn_ass_02(tr_as: &Vec<PeaAssVar>, fnm: &str) -> Result<String, Box<dyn Error>> {
    let flds = [
        VarType::NewCarRegVp01,
        VarType::GppVp02,
        VarType::MaxPosPowSubVs01,
        VarType::MaxNegPowSubVs02,
        VarType::VsppMvVs03,
        VarType::SppHvVs04,
        VarType::BigLotMvVs05,
        VarType::BigLotHvVs06,
        VarType::SubPowCapVs07,
        VarType::MaxPosPowFeederVf01,
        VarType::MaxNegPowFeederVf02,
        VarType::MaxPosDiffFeederVf03,
        VarType::MaxNegDiffFeederVf04,
        VarType::NoMeterTransVt01,
        VarType::SmallSellTrVt02,
        VarType::ChgStnCapTrVt03,
        VarType::ChgStnSellTrVt04,
        VarType::PwCapTriVt05,
        VarType::ZoneTrVt06,
        VarType::PopTrVt07,
        VarType::UnbalPowTrVt08,
        VarType::PkPowTrVt09,
        VarType::LargeSellTrVt10,
        VarType::AllNoMeterTrVt11,
        VarType::NoTrVt12,
        VarType::HmChgEvTrVc01,
        VarType::LvPowSatTrVc02,
        VarType::CntLvPowSatTrVc03,
        VarType::ChgStnCapVc04,
        VarType::ChgStnSellVc05,
        VarType::MvPowSatTrVc06,
        VarType::PowSolarVc07,
        VarType::MvVsppVc08,
        VarType::HvSppVc09,
        VarType::SmallSellVc10,
        VarType::LargeSellVc11,
        VarType::UnbalPowVc12,
        VarType::CntUnbalPowVc13,
        VarType::Uc1ValVc14,
        VarType::Uc2ValVc15,
        VarType::Uc3ValVc16,
        VarType::Uc1RankVc17,
        VarType::Uc2RankVc18,
        VarType::Uc3RankVc19,
        VarType::NoHmChgEvTr,
        VarType::PowHmChgEvTr,
        VarType::PkSelPowPhsAKw,
        VarType::PkSelPowPhsBKw,
        VarType::PkSelPowPhsCKw,
        VarType::PkSelPowPhsAvg,
        VarType::PkSelPowPhsMax,
        VarType::UnbalPowRate,
        VarType::TransLossKw,
        VarType::UnbalPowLossKw,
        VarType::TakeNote,
        VarType::CntTrUnbalLoss,
        VarType::CntTrSatLoss,
        VarType::SolarEnergy,
    ];
    write_text_02(tr_as, &flds, fnm)
}

fn write_text_02(
    tr_as: &Vec<PeaAssVar>,
    flds: &[VarType],
    fnm: &str,
) -> Result<String, Box<dyn Error>> {
    let mut x = String::new();
    use std::fmt::Write;
    for t in tr_as {
        //t._a("y:");
        write!(x, "{}", t.sbid)?;
        write!(x, "\t{}", t.pvid)?;
        write!(x, "\t{}", t.arid)?;
        for f in flds.iter() {
            let d = t.v[f.clone() as usize].v;
            write!(x, "\t{d}")?;
        }
        writeln!(x)?;
    }
    println!("        ===== write to {fnm}");
    let b = x.as_bytes();
    let h = sha256::digest(b);
    std::fs::write(fnm, b)?;
    Ok(h)
}

fn write_text_01(
    tr_as: &Vec<PeaAssVar>,
    flds: &[VarType],
    fnm: &str,
) -> Result<String, Box<dyn Error>> {
    let mut x = String::new();
    use std::fmt::Write;
    for t in tr_as {
        //t._a("y:");
        write!(x, "{}", t.sbid)?;
        write!(x, "\t{}", t.fdid)?;
        write!(x, "\t{}", t.aoj)?;
        write!(x, "\t{}", t.own)?;
        write!(x, "\t{}", t.peano)?;
        for f in flds.iter() {
            let d = t.v[f.clone() as usize].v;
            write!(x, "\t{d}")?;
        }
        writeln!(x)?;
    }
    println!("        ===== write to {fnm}");
    let b = x.as_bytes();
    let h = sha256::digest(b);
    std::fs::write(fnm, b)?;
    Ok(h)
}

/// ประมวลผลรวมเพื่อเกณฑ์การคัดเลือก
/// summery transformaters to substation
pub fn c01_chk_03() -> Result<(), Box<dyn Error>> {
    let dnm = "/mnt/e/CHMBACK/pea-data/c01_pea";
    let buf = std::fs::read(format!("{dnm}/000_pea.bin")).unwrap();
    let (pea, _): (Pea, usize) =
        bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
    let mut aids: Vec<_> = pea.aream.keys().collect();
    aids.sort();
    let subhs = p01_chk();
    //
    //let mut pvcn = 0;
    let mut v_pvas = Vec::<PeaAssVar>::new();
    let mut v_sbas = Vec::<PeaAssVar>::new();
    let mut sbas_mx = PeaAssVar::default();
    for aid in aids {
        let Some(ar) = pea.aream.get(aid) else {
            continue;
        };
        let mut pids: Vec<_> = ar.provm.keys().collect();
        pids.sort();
        for pid in pids {
            let Some(prov) = ar.provm.get(pid) else {
                continue;
            };
            let mut pvas = PeaAssVar::from(0u64);
            pvas.arid = aid.to_string();
            pvas.pvid = pid.to_string();
            //println!("  pv:{pid}");
            let mut sids: Vec<_> = prov.subm.keys().collect();
            sids.sort();
            for sid in sids {
                let Some(_sb) = prov.subm.get(sid) else {
                    continue;
                };
                let Ok(buf) = std::fs::read(format!("{dnm}/{sid}-rw3.bin")) else {
                    continue;
                };
                let (v_tras_raw, _): (Vec<PeaAssVar>, usize) =
                    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
                if v_tras_raw.is_empty() {
                    println!("    {sid} - NO data ");
                    continue;
                }
                let tras = &v_tras_raw[0];
                let mut sbas = PeaAssVar::from(0u64);
                sbas.arid = aid.to_string();
                sbas.pvid = pid.to_string();
                sbas.sbid = tras.sbid.to_string();
                let note = if subhs.contains(&sbas.sbid) {
                    1f32
                } else {
                    0f32
                };

                for tras in &v_tras_raw {
                    sbas.add(tras);
                }
                sbas.copy(tras, VarType::NewCarRegVp01);
                sbas.copy(tras, VarType::GppVp02);
                sbas.copy(tras, VarType::MaxPosPowSubVs01);
                sbas.copy(tras, VarType::MaxNegPowSubVs02);
                sbas.copy(tras, VarType::VsppMvVs03);
                sbas.copy(tras, VarType::SppHvVs04);
                sbas.copy(tras, VarType::BigLotMvVs05);
                sbas.copy(tras, VarType::BigLotHvVs06);
                sbas.copy(tras, VarType::SubPowCapVs07);
                sbas.copy(tras, VarType::SolarEnergy);
                let solar = sbas.v[VarType::SolarEnergy as usize].v;
                if solar > 0f32 {
                    println!(">>>>>>>>>>> {sid} solar: {solar} =============");
                }

                sbas.v[VarType::TakeNote as usize].v = note;
                sbas_mx.max(&sbas);
                pvas.add(&sbas);
                pvas.copy(tras, VarType::NewCarRegVp01);
                pvas.copy(tras, VarType::GppVp02);
                v_sbas.push(sbas);
                //println!("   {sid} - {}", v_tras.len());
            } // end sub loop
            v_pvas.push(pvas);
        } // end provi loop
    } // end area
    let mut uc1_v: Vec<_> = v_sbas
        .iter()
        .enumerate()
        .map(|(i, s)| (s.v[VarType::Uc1ValVc14 as usize].v, i))
        .collect();
    uc1_v.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    for (r, (_, i)) in uc1_v.iter().enumerate() {
        v_sbas[*i].v[VarType::Uc1RankVc17 as usize].v = r as f32;
    }

    let mut uc2_v: Vec<_> = v_sbas
        .iter()
        .enumerate()
        .map(|(i, s)| (s.v[VarType::Uc2ValVc15 as usize].v, i))
        .collect();
    uc2_v.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    for (r, (_, i)) in uc2_v.iter().enumerate() {
        v_sbas[*i].v[VarType::Uc2RankVc18 as usize].v = r as f32;
    }

    let mut uc3_v: Vec<_> = v_sbas
        .iter()
        .enumerate()
        .map(|(i, s)| (s.v[VarType::Uc3ValVc16 as usize].v, i))
        .collect();
    uc3_v.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    for (r, (_, i)) in uc3_v.iter().enumerate() {
        v_sbas[*i].v[VarType::Uc3RankVc19 as usize].v = r as f32;
    }

    // save ev bin
    let bin: Vec<u8> = bincode::encode_to_vec(&v_sbas, bincode::config::standard()).unwrap();
    std::fs::write(format!("{dnm}/000-sbrw.bin"), bin).unwrap();
    write_trn_ass_02(&v_sbas, &format!("{dnm}/000-sbrw0.txt"))?;
    write_ass_csv_02(&v_sbas, &format!("{dnm}/000-sbrw0.csv"))?;

    let mut v_sbas_no = v_sbas.clone();
    for sub in v_sbas_no.iter_mut() {
        sub.nor(&sbas_mx);
    }
    let bin: Vec<u8> = bincode::encode_to_vec(&v_sbas_no, bincode::config::standard()).unwrap();
    std::fs::write(format!("{dnm}/000-sbno.bin"), bin).unwrap();
    write_trn_ass_02(&v_sbas_no, &format!("{dnm}/000-sbno0.txt"))?;
    write_ass_csv_02(&v_sbas_no, &format!("{dnm}/000-sbno0.csv"))?;

    let bin: Vec<u8> = bincode::encode_to_vec(&v_pvas, bincode::config::standard()).unwrap();
    std::fs::write(format!("{dnm}/000-pvrw.bin"), bin).unwrap();
    write_trn_ass_02(&v_pvas, &format!("{dnm}/000-pvrw.txt"))?;
    write_ass_csv_02(&v_sbas_no, &format!("{dnm}/000-pvrw.csv"))?;
    Ok(())
}

/// check for ราชบุรี if it belongs to multiple area
pub fn c01_chk_04() -> Result<(), Box<dyn Error>> {
    let g0 = ProcEngine::prep1();
    let mut pea = Pea::default();
    let mut pvs = HashSet::<String>::new();
    for (sb, sf) in &g0.sbif {
        let ar = sf.arid.to_string();
        let ar_e = pea.aream.entry(ar).or_insert_with(|| PeaArea {
            arid: sf.arid.to_string(),
            ..Default::default()
        });
        let pv_e = ar_e
            .provm
            .entry(sf.prov.to_string())
            .or_insert_with(|| PeaProv {
                pvnm: sf.prov.to_string(),
                ..Default::default()
            });
        let _sb_e = pv_e.subm.entry(sb.to_string()).or_insert_with(|| PeaSub {
            sbid: sb.to_string(),
            name: sf.name.clone(),
            enam: sf.enam.clone(),
            area: sf.area.clone(),
            arid: sf.arid.clone(),
            prov: sf.prov.clone(),
            ..Default::default()
        });

        pvs.insert(sf.prov.to_string());
    }
    for (_aid, ar) in &pea.aream {
        //println!("{}", ar.arid);
        for (pv, pf) in &ar.provm {
            if pv != "ราชบุรี" {
                continue;
            }
            println!("{pv}");
            for (sid, sb) in &pf.subm {
                println!(
                    "  {sid} n:{} e:{} a:{}-{} p:{}",
                    sb.name, sb.enam, sb.area, sb.arid, sb.prov
                );
            }
        }
    }

    Ok(())
}

fn write_ass_csv_01(tr_as: &Vec<PeaAssVar>, fnm: &str) -> Result<String, Box<dyn Error>> {
    let mut x = String::new();
    use std::fmt::Write;
    let flds = [
        VarType::NewCarRegVp01,
        VarType::GppVp02,
        VarType::MaxPosPowSubVs01,
        VarType::MaxNegPowSubVs02,
        VarType::VsppMvVs03,
        VarType::SppHvVs04,
        VarType::BigLotMvVs05,
        VarType::BigLotHvVs06,
        VarType::SubPowCapVs07,
        VarType::MaxPosPowFeederVf01,
        VarType::MaxNegPowFeederVf02,
        VarType::MaxPosDiffFeederVf03,
        VarType::MaxNegDiffFeederVf04,
        VarType::NoMeterTransVt01,
        VarType::SmallSellTrVt02,
        VarType::ChgStnCapTrVt03,
        VarType::ChgStnSellTrVt04,
        VarType::PwCapTriVt05,
        VarType::ZoneTrVt06,
        VarType::PopTrVt07,
        VarType::UnbalPowTrVt08,
        VarType::PkPowTrVt09,
        VarType::LargeSellTrVt10,
        VarType::AllNoMeterTrVt11,
        VarType::NoTrVt12,
        VarType::HmChgEvTrVc01,
        VarType::LvPowSatTrVc02,
        VarType::CntLvPowSatTrVc03,
        VarType::ChgStnCapVc04,
        VarType::ChgStnSellVc05,
        VarType::MvPowSatTrVc06,
        VarType::PowSolarVc07,
        VarType::MvVsppVc08,
        VarType::HvSppVc09,
        VarType::SmallSellVc10,
        VarType::LargeSellVc11,
        VarType::UnbalPowVc12,
        VarType::CntUnbalPowVc13,
        VarType::Uc1ValVc14,
        VarType::Uc2ValVc15,
        VarType::Uc3ValVc16,
        VarType::Uc1RankVc17,
        VarType::Uc2RankVc18,
        VarType::Uc3RankVc19,
        VarType::NoHmChgEvTr,
        VarType::PowHmChgEvTr,
    ];
    write!(x, "\"SUB\"")?;
    write!(x, ",\"FEEDER\"")?;
    write!(x, ",\"AOJ\"")?;
    write!(x, ",\"OWN\"")?;
    write!(x, ",\"PEANO\"")?;
    for f in flds.iter() {
        let l = format!("{f:?}");
        write!(x, ",\"{l}\"")?;
    }
    writeln!(x)?;
    for t in tr_as {
        //t._a("y:");
        write!(x, "\"{}\"", t.sbid)?;
        write!(x, ",\"{}\"", t.fdid)?;
        write!(x, ",\"{}\"", t.aoj)?;
        write!(x, ",\"{}\"", t.own)?;
        write!(x, ",\"{}\"", t.peano)?;
        write!(x, ",\"{}\"", t.v[VarType::NewCarRegVp01 as usize].v)?;
        for f in flds.iter() {
            let d = t.v[f.clone() as usize].v;
            write!(x, ",{d}")?;
        }
        writeln!(x)?;
    }
    println!("        ===== write to {fnm}");
    let b = x.as_bytes();
    let h = sha256::digest(b);
    std::fs::write(fnm, b)?;

    Ok(h)
}

fn write_ass_csv_02(tr_as: &Vec<PeaAssVar>, fnm: &str) -> Result<String, Box<dyn Error>> {
    let mut x = String::new();
    use std::fmt::Write;
    let flds = [
        VarType::NewCarRegVp01,
        VarType::GppVp02,
        VarType::MaxPosPowSubVs01,
        VarType::MaxNegPowSubVs02,
        VarType::VsppMvVs03,
        VarType::SppHvVs04,
        VarType::BigLotMvVs05,
        VarType::BigLotHvVs06,
        VarType::SubPowCapVs07,
        VarType::MaxPosPowFeederVf01,
        VarType::MaxNegPowFeederVf02,
        VarType::MaxPosDiffFeederVf03,
        VarType::MaxNegDiffFeederVf04,
        VarType::NoMeterTransVt01,
        VarType::SmallSellTrVt02,
        VarType::ChgStnCapTrVt03,
        VarType::ChgStnSellTrVt04,
        VarType::PwCapTriVt05,
        VarType::ZoneTrVt06,
        VarType::PopTrVt07,
        VarType::UnbalPowTrVt08,
        VarType::PkPowTrVt09,
        VarType::LargeSellTrVt10,
        VarType::AllNoMeterTrVt11,
        VarType::NoTrVt12,
        VarType::HmChgEvTrVc01,
        VarType::LvPowSatTrVc02,
        VarType::CntLvPowSatTrVc03,
        VarType::ChgStnCapVc04,
        VarType::ChgStnSellVc05,
        VarType::MvPowSatTrVc06,
        VarType::PowSolarVc07,
        VarType::MvVsppVc08,
        VarType::HvSppVc09,
        VarType::SmallSellVc10,
        VarType::LargeSellVc11,
        VarType::UnbalPowVc12,
        VarType::CntUnbalPowVc13,
        VarType::Uc1ValVc14,
        VarType::Uc2ValVc15,
        VarType::Uc3ValVc16,
        VarType::Uc1RankVc17,
        VarType::Uc2RankVc18,
        VarType::Uc3RankVc19,
        VarType::NoHmChgEvTr,
        VarType::PowHmChgEvTr,
        VarType::SolarEnergy,
    ];
    write!(x, "\"SUB\"")?;
    write!(x, ",\"PROV\"")?;
    write!(x, ",\"ARID\"")?;
    for f in flds.iter() {
        let l = format!("{f:?}");
        write!(x, ",\"{l}\"")?;
    }
    writeln!(x)?;
    for t in tr_as {
        //t._a("y:");
        write!(x, "\"{}\"", t.sbid)?;
        write!(x, ",\"{}\"", t.pvid)?;
        write!(x, ",\"{}\"", t.arid)?;
        for f in flds.iter() {
            let d = t.v[f.clone() as usize].v;
            write!(x, ",{d}")?;
        }
        writeln!(x)?;
    }
    println!("        ===== write to {fnm}");
    let b = x.as_bytes();
    let h = sha256::digest(b);
    std::fs::write(fnm, b)?;

    Ok(h)
}

pub static EV_LIKELY: phf::Map<&'static str, f32> = phf_map! {
"ระยอง" => 1f32,
"ชลบุรี" => 1f32,
"ปทุมธานี" => 1f32,
"สมุทรสาคร" => 1f32,
"นครปฐม" => 1f32,
"สงขลา" => 1f32,
"พระนครศรีอยุธยา" => 1f32,
"สระบุรี" => 1f32,
"เชียงใหม่" => 1f32,
"ฉะเชิงเทรา" => 1f32,
"นครราชสีมา" => 1f32,
"ราชบุรี" => 1f32,
"ขอนแก่น" => 1f32,
"ปราจีนบุรี" => 1f32,
"พิษณุโลก" => 1f32,
"สุราษฎร์ธานี" => 1f32,
"นครสวรรค์" => 1f32,
"เพชรบุรี" => 1f32,
"ภูเก็ต" => 1f32,
};

pub static SELE_LIKELY: phf::Map<&'static str, f32> = phf_map! {
"ระยอง" => 1f32,
"ชลบุรี" => 1f32,
"ปทุมธานี" => 1f32,
"สมุทรสาคร" => 1f32,
"นครปฐม" => 1f32,
"สงขลา" => 1f32,
"พระนครศรีอยุธยา" => 1f32,
"สระบุรี" => 1f32,
"เชียงใหม่" => 1f32,
"ฉะเชิงเทรา" => 1f32,
"นครราชสีมา" => 1f32,
"ราชบุรี" => 1f32,
"ขอนแก่น" => 1f32,
"ปราจีนบุรี" => 1f32,
"พิษณุโลก" => 1f32,
"สุราษฎร์ธานี" => 1f32,
"นครสวรรค์" => 1f32,
"เพชรบุรี" => 1f32,
"ภูเก็ต" => 1f32,
};

/*
"สมุทรปราการ" => 320294,
"กระบี่" => 174058,
"สมุทรสงคราม" => 167164,
"ลพบุรี" => 152831,
"สระแก้ว" => 82526,
"บุรีรัมย์" => 91636,
"ระนอง" => 99331,
"จันทบุรี" => 253522,
"ลำพูน" => 236619,
"ชุมพร" => 230319,
"พังงา" => 229213,
"ประจวบคีรีขันธ์" => 221151,
"นนทบุรี" => 214515,
"ตราด" => 164835,
"ชัยนาท" => 157159,
"กำแพงเพชร" => 155404,
"กาญจนบุรี" => 153662,
"สิงห์บุรี" => 151750,
"อ่างทอง" => 135248,
"นครศรีธรรมราช" => 127405,
"นครนายก" => 126435,
"สุพรรณบุรี" => 124482,
"อุทัยธานี" => 123946,
"ตาก" => 121537,
"อุตรดิตถ์" => 120720,
"เลย" => 117624,
"ตรัง" => 111746,
"สตูล" => 110312,
"พะเยา" => 109275,
"ยะลา" => 108108,
"ลำปาง" => 107732,
"หนองคาย" => 107589,
"พิจิตร" => 105054,
"เชียงราย" => 102988,
"เพชรบูรณ์" => 100936,
"อุดรธานี" => 100005,
"นครพนม" => 96731,
"สุโขทัย" => 93208,
"แพร่" => 91324,
"ศรีสะเกษ" => 91060,
"มหาสารคาม" => 90996,
"สุรินทร์" => 89852,
"น่าน" => 89515,
"พัทลุง" => 87098,
"ชัยภูมิ" => 85951,
"อำนาจเจริญ" => 85707,
"กาฬสินธุ์" => 84785,
"บึงกาฬ" => 84021,
"ปัตตานี" => 83369,
"อุบลราชธานี" => 82895,
"ร้อยเอ็ด" => 82491,
"สกลนคร" => 78895,
"ยโสธร" => 77376,
"มุกดาหาร" => 72251,
"แม่ฮ่องสอน" => 69828,
"หนองบัวลำภู" => 69008,
"นราธิวาส" => 64005,
*/
