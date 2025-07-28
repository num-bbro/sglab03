/*
use crate::p01::get_tr_den;
use crate::p01::get_tr_sorf;
use crate::p01::get_tr_volta;
use crate::p01::get_tr_zone;
use crate::p01::mon_kwh_2_kw;
use crate::p01::p01_chk;
use crate::p01::trf_kva_2_kw;
use crate::p01::AojObj;
use sglab02_lib::sg::gis1::ar_list;
use sglib04::geo4::GPPS;
use crate::p01::SubAssObj;
*/
use crate::p01::ProcEngine;
use bincode::{Decode, Encode};
use regex::Regex;
use std::collections::HashMap;
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
pub struct PeaAssVar {
    pub arid: String,
    pub pvid: String,
    pub sbid: String,
    pub fdid: String,
    pub n1d: String,
    pub own: String,
    pub peano: String,
    pub aoj: String,
    pub set: bool,
    //  province level
    pub vp01: f32, // EV ratio of province
    pub vp02: f32, // GPP
    // substation level
    pub vs01: f32,
    pub vs02: f32,
    pub vs03: f32,
    pub vs04: f32,
    pub vs05: f32,
    pub vs06: f32,
    pub vs07: f32,
    // feeder level
    pub vf01: f32,
    pub vf02: f32,
    pub vf03: f32,
    pub vf04: f32,
    // transformer level
    pub vt01: f32, // meter
    pub vt02: f32, // sold energy May24
    pub vt03: f32, // charge station : kw
    pub vt04: f32,
    pub vt05: f32, // transformer capacity : kw
    pub vt06: f32, // zone 1-5
    pub vt07: f32, // population density
    pub vt08: f32,
    pub vt09: f32,
    pub vt10: f32,
    pub vt11: f32,
    pub vt12: f32,
    // calculate
    pub vc01: f32,
    pub vc02: f32,
    pub vc03: f32,
    pub vc04: f32,
    pub vc05: f32,
    pub vc06: f32,
    pub vc07: f32,
    pub vc08: f32,
    pub vc09: f32,
    pub vc10: f32,
    pub vc11: f32,
    pub vc12: f32,
    // subm
    pub sum: f32,
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
        self.vp01 /= o;
        self.vp02 /= o;
        self.vs01 /= o;
        self.vs02 /= o;
        self.vs03 /= o;
        self.vs04 /= o;
        self.vs05 /= o;
        self.vs06 /= o;
        self.vs07 /= o;
        self.vf01 /= o;
        self.vf02 /= o;
        self.vf03 /= o;
        self.vf04 /= o;
        self.vt01 /= o;
        self.vt02 /= o;
        self.vt03 /= o;
        self.vt04 /= o;
        self.vt05 /= o;
        self.vt06 /= o;
        self.vt07 /= o;
        self.vt08 /= o;
        self.vt09 /= o;
        self.vt10 /= o;
        self.vt11 /= o;
        self.vt12 /= o;
        self.vc01 /= o;
        self.vc02 /= o;
        self.vc03 /= o;
        self.vc04 /= o;
        self.vc05 /= o;
        self.vc06 /= o;
        self.vc07 /= o;
        self.vc08 /= o;
        self.vc09 /= o;
        self.vc10 /= o;
        self.vc11 /= o;
        self.vc12 /= o;
    }
    pub fn nor(&mut self, o: &PeaAssVar) {
        self.vp01 /= z2o(o.vp01);
        self.vp02 /= z2o(o.vp02);
        self.vs01 /= z2o(o.vs01);
        self.vs02 /= z2o(o.vs02);
        self.vs03 /= z2o(o.vs03);
        self.vs04 /= z2o(o.vs04);
        self.vs05 /= z2o(o.vs05);
        self.vs06 /= z2o(o.vs06);
        self.vs07 /= z2o(o.vs07);
        self.vf01 /= z2o(o.vf01);
        self.vf02 /= z2o(o.vf02);
        self.vf03 /= z2o(o.vf03);
        self.vf04 /= z2o(o.vf04);
        self.vt01 /= z2o(o.vt01);
        self.vt02 /= z2o(o.vt02);
        self.vt03 /= z2o(o.vt03);
        self.vt04 /= z2o(o.vt04);
        self.vt05 /= z2o(o.vt05);
        self.vt06 /= z2o(o.vt06);
        self.vt07 /= z2o(o.vt07);
        self.vt08 /= z2o(o.vt08);
        self.vt09 /= z2o(o.vt09);
        self.vt10 /= z2o(o.vt10);
        self.vt11 /= z2o(o.vt11);
        self.vt12 /= z2o(o.vt12);
        self.vc01 /= z2o(o.vc01);
        self.vc02 /= z2o(o.vc02);
        self.vc03 /= z2o(o.vc03);
        self.vc04 /= z2o(o.vc04);
        self.vc05 /= z2o(o.vc05);
        self.vc06 /= z2o(o.vc06);
        self.vc07 /= z2o(o.vc07);
        self.vc08 /= z2o(o.vc08);
        self.vc09 /= z2o(o.vc09);
        self.vc10 /= z2o(o.vc10);
        self.vc11 /= z2o(o.vc11);
        self.vc12 /= z2o(o.vc12);
    }
    pub fn add(&mut self, o: &PeaAssVar) {
        self.vp01 += o.vp01;
        self.vp02 += o.vp02;
        self.vs01 += o.vs01;
        self.vs02 += o.vs02;
        self.vs03 += o.vs03;
        self.vs04 += o.vs04;
        self.vs05 += o.vs05;
        self.vs06 += o.vs06;
        self.vs07 += o.vs07;
        self.vf01 += o.vf01;
        self.vf02 += o.vf02;
        self.vf03 += o.vf03;
        self.vf04 += o.vf04;
        self.vt01 += o.vt01;
        self.vt02 += o.vt02;
        self.vt03 += o.vt03;
        self.vt04 += o.vt04;
        self.vt05 += o.vt05;
        self.vt06 += o.vt06;
        self.vt07 += o.vt07;
        self.vt08 += o.vt08;
        self.vt09 += o.vt09;
        self.vt10 += o.vt10;
        self.vt11 += o.vt11;
        self.vt12 += o.vt12;
        self.vc01 += o.vc01;
        self.vc02 += o.vc02;
        self.vc03 += o.vc03;
        self.vc04 += o.vc04;
        self.vc05 += o.vc05;
        self.vc06 += o.vc06;
        self.vc07 += o.vc07;
        self.vc08 += o.vc08;
        self.vc09 += o.vc09;
        self.vc10 += o.vc10;
        self.vc11 += o.vc11;
        self.vc12 += o.vc12;
    }
    pub fn max(&mut self, o: &PeaAssVar) {
        if !self.set {
            self.set = true;
            self.vp01 = o.vp01;
            self.vp02 = o.vp02;
            self.vs01 = o.vs01;
            self.vs02 = o.vs02;
            self.vs03 = o.vs03;
            self.vs04 = o.vs04;
            self.vs05 = o.vs05;
            self.vs06 = o.vs06;
            self.vs07 = o.vs07;
            self.vf01 = o.vf01;
            self.vf02 = o.vf02;
            self.vf03 = o.vf03;
            self.vf04 = o.vf04;
            self.vt01 = o.vt01;
            self.vt02 = o.vt02;
            self.vt03 = o.vt03;
            self.vt04 = o.vt04;
            self.vt05 = o.vt05;
            self.vt06 = o.vt06;
            self.vt07 = o.vt07;
            self.vt08 = o.vt08;
            self.vt09 = o.vt09;
            self.vt10 = o.vt10;
            self.vt11 = o.vt11;
            self.vt12 = o.vt12;
            self.vc01 = o.vc01;
            self.vc02 = o.vc02;
            self.vc03 = o.vc03;
            self.vc04 = o.vc04;
            self.vc05 = o.vc05;
            self.vc06 = o.vc06;
            self.vc07 = o.vc07;
            self.vc08 = o.vc08;
            self.vc09 = o.vc09;
            self.vc10 = o.vc10;
            self.vc11 = o.vc11;
            self.vc12 = o.vc12;
            return;
        }
        self.vp01 = self.vp01.max(o.vp01);
        self.vp02 = self.vp02.max(o.vp02);
        self.vs01 = self.vs01.max(o.vs01);
        self.vs02 = self.vs02.max(o.vs02);
        self.vs03 = self.vs03.max(o.vs03);
        self.vs04 = self.vs04.max(o.vs04);
        self.vs05 = self.vs05.max(o.vs05);
        self.vs06 = self.vs06.max(o.vs06);
        self.vs07 = self.vs07.max(o.vs07);
        self.vf01 = self.vf01.max(o.vf01);
        self.vf02 = self.vf02.max(o.vf02);
        self.vf03 = self.vf03.max(o.vf03);
        self.vf04 = self.vf04.max(o.vf04);
        self.vt01 = self.vt01.max(o.vt01);
        self.vt02 = self.vt02.max(o.vt02);
        self.vt03 = self.vt03.max(o.vt03);
        self.vt04 = self.vt04.max(o.vt04);
        self.vt05 = self.vt05.max(o.vt05);
        self.vt06 = self.vt06.max(o.vt06);
        self.vt07 = self.vt07.max(o.vt07);
        self.vt08 = self.vt08.max(o.vt08);
        self.vt09 = self.vt09.max(o.vt09);
        self.vt10 = self.vt10.max(o.vt10);
        self.vt11 = self.vt11.max(o.vt11);
        self.vt12 = self.vt12.max(o.vt12);
        self.vc01 = self.vc01.max(o.vc01);
        self.vc02 = self.vc02.max(o.vc02);
        self.vc03 = self.vc03.max(o.vc03);
        self.vc04 = self.vc04.max(o.vc04);
        self.vc05 = self.vc05.max(o.vc05);
        self.vc06 = self.vc06.max(o.vc06);
        self.vc07 = self.vc07.max(o.vc07);
        self.vc08 = self.vc08.max(o.vc08);
        self.vc09 = self.vc09.max(o.vc09);
        self.vc10 = self.vc10.max(o.vc10);
        self.vc11 = self.vc11.max(o.vc11);
        self.vc12 = self.vc12.max(o.vc12);
    }
    pub fn min(&mut self, o: &PeaAssVar) {
        if !self.set {
            self.set = true;
            self.vp01 = o.vp01;
            self.vp02 = o.vp02;
            self.vs01 = o.vs01;
            self.vs02 = o.vs02;
            self.vs03 = o.vs03;
            self.vs04 = o.vs04;
            self.vs05 = o.vs05;
            self.vs06 = o.vs06;
            self.vs07 = o.vs07;
            self.vf01 = o.vf01;
            self.vf02 = o.vf02;
            self.vf03 = o.vf03;
            self.vf04 = o.vf04;
            self.vt01 = o.vt01;
            self.vt02 = o.vt02;
            self.vt03 = o.vt03;
            self.vt04 = o.vt04;
            self.vt05 = o.vt05;
            self.vt06 = o.vt06;
            self.vt07 = o.vt07;
            self.vt08 = o.vt08;
            self.vt09 = o.vt09;
            self.vt10 = o.vt10;
            self.vt11 = o.vt11;
            self.vt12 = o.vt12;
            self.vc01 = o.vc01;
            self.vc02 = o.vc02;
            self.vc03 = o.vc03;
            self.vc04 = o.vc04;
            self.vc05 = o.vc05;
            self.vc06 = o.vc06;
            self.vc07 = o.vc07;
            self.vc08 = o.vc08;
            self.vc09 = o.vc09;
            self.vc10 = o.vc10;
            self.vc11 = o.vc11;
            self.vc12 = o.vc12;
            return;
        }
        self.vp01 = self.vp01.min(o.vp01);
        self.vp02 = self.vp02.min(o.vp02);
        self.vs01 = self.vs01.min(o.vs01);
        self.vs02 = self.vs02.min(o.vs02);
        self.vs03 = self.vs03.min(o.vs03);
        self.vs04 = self.vs04.min(o.vs04);
        self.vs05 = self.vs05.min(o.vs05);
        self.vs06 = self.vs06.min(o.vs06);
        self.vs07 = self.vs07.min(o.vs07);
        self.vf01 = self.vf01.min(o.vf01);
        self.vf02 = self.vf02.min(o.vf02);
        self.vf03 = self.vf03.min(o.vf03);
        self.vf04 = self.vf04.min(o.vf04);
        self.vt01 = self.vt01.min(o.vt01);
        self.vt02 = self.vt02.min(o.vt02);
        self.vt03 = self.vt03.min(o.vt03);
        self.vt04 = self.vt04.min(o.vt04);
        self.vt05 = self.vt05.min(o.vt05);
        self.vt06 = self.vt06.min(o.vt06);
        self.vt07 = self.vt07.min(o.vt07);
        self.vt08 = self.vt08.min(o.vt08);
        self.vt09 = self.vt09.min(o.vt09);
        self.vt10 = self.vt10.min(o.vt10);
        self.vt11 = self.vt11.min(o.vt11);
        self.vt12 = self.vt12.min(o.vt12);
        self.vc01 = self.vc01.min(o.vc01);
        self.vc02 = self.vc02.min(o.vc02);
        self.vc03 = self.vc03.min(o.vc03);
        self.vc04 = self.vc04.min(o.vc04);
        self.vc05 = self.vc05.min(o.vc05);
        self.vc06 = self.vc06.min(o.vc06);
        self.vc07 = self.vc07.min(o.vc07);
        self.vc08 = self.vc08.min(o.vc08);
        self.vc09 = self.vc09.min(o.vc09);
        self.vc10 = self.vc10.min(o.vc10);
        self.vc11 = self.vc11.min(o.vc11);
        self.vc12 = self.vc12.min(o.vc12);
    }
    pub fn weigh(&mut self, o: &PeaAssVar) {
        self.vp01 *= o.vp01;
        self.vp02 *= o.vp02;
        self.vs01 *= o.vs01;
        self.vs02 *= o.vs02;
        self.vs03 *= o.vs03;
        self.vs04 *= o.vs04;
        self.vs05 *= o.vs05;
        self.vs06 *= o.vs06;
        self.vs07 *= o.vs07;
        self.vf01 *= o.vf01;
        self.vf02 *= o.vf02;
        self.vf03 *= o.vf03;
        self.vf04 *= o.vf04;
        self.vt01 *= o.vt01;
        self.vt02 *= o.vt02;
        self.vt03 *= o.vt03;
        self.vt04 *= o.vt04;
        self.vt05 *= o.vt05;
        self.vt06 *= o.vt06;
        self.vt07 *= o.vt07;
        self.vt08 *= o.vt08;
        self.vt09 *= o.vt09;
        self.vt10 *= o.vt10;
        self.vt11 *= o.vt11;
        self.vt12 *= o.vt12;
        self.vc01 *= o.vc01;
        self.vc02 *= o.vc02;
        self.vc03 *= o.vc03;
        self.vc04 *= o.vc04;
        self.vc05 *= o.vc05;
        self.vc06 *= o.vc06;
        self.vc07 *= o.vc07;
        self.vc08 *= o.vc08;
        self.vc09 *= o.vc09;
        self.vc10 *= o.vc10;
        self.vc11 *= o.vc11;
        self.vc12 *= o.vc12;
    }
    pub fn sum(&mut self) {
        self.sum = self.vp01
            + self.vp02
            + self.vs01
            + self.vs02
            + self.vs03
            + self.vs04
            + self.vs05
            + self.vs06
            + self.vs07
            + self.vf01
            + self.vf02
            + self.vf03
            + self.vf04
            + self.vt01
            + self.vt02
            + self.vt03
            + self.vt04
            + self.vt05
            + self.vt06
            + self.vt07
            + self.vt08
            + self.vt09
            + self.vt10
            + self.vt11
            + self.vt12
            + self.vc01
            + self.vc02
            + self.vc03
            + self.vc04
            + self.vc05
            + self.vc06
            + self.vc07
            + self.vc08
            + self.vc09
            + self.vc10
            + self.vc11
            + self.vc12;
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
    std::fs::create_dir_all(dnm)?;
    let smrt = Regex::new(r"[12].*").unwrap();
    let g0 = ProcEngine::prep1();
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
                    println!("     ================= NO1 {sid}");
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
                    println!("     ================= NO2 {sid}");
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
                            let tp = if smrt.captures(bl.rate.as_str()).is_some()
                                && bl.main.is_empty()
                            {
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
                    let Some(ref code) = ao.code else {
                        continue;
                    };
                    let Some(ref sht_name) = ao.sht_name else {
                        continue;
                    };
                    let Some(ref office) = ao.office else {
                        continue;
                    };
                    let Some(ref pea) = ao.pea else {
                        continue;
                    };
                    let Some(ref aoj_sz) = ao.aoj_sz else {
                        continue;
                    };
                    let Some(ref reg) = ao.reg else {
                        continue;
                    };
                    let Some(ref name) = ao.name else {
                        continue;
                    };
                    let Some(ref level) = ao.level else {
                        continue;
                    };
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
                        level: *level,
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
    let mut cn0 = 0;
    let mut cn1 = 1;
    let (mut mtc0, mut mtz0) = (0, 0);
    let (mut txpc, mut txcc) = (0, 0);

    let evwe = PeaAssVar {
        vp01: 0.20,
        vp02: 0.05,
        vt01: 0.20,
        vt02: 0.20,
        vt05: 0.20,
        vt06: 0.05,
        vt07: 0.05,
        ..Default::default()
    };
    let sowe = PeaAssVar {
        vp02: 0.10,
        vs01: 0.02,
        vs02: 0.02,
        vs03: 0.02,
        vs04: 0.02,
        vs05: 0.02,
        vf02: 0.05,
        vf03: 0.05,
        vt01: 0.20,
        vt02: 0.20,
        vt05: 0.20,
        vt06: 0.05,
        vt07: 0.05,
        ..Default::default()
    };

    let mut tras_mx1 = PeaAssVar::default();
    let e0 = ProcEngine::prep1();

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
    for id in &aids {
        let aid = id.to_string();
        let Some(ar) = pea.aream.get(&aid) else {
            continue;
        };
        println!("ar:{}", ar.arid);
        //let eg = ProcEngine::prep3(id);
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

            println!("  pv:{pid}");
            let mut sids: Vec<_> = prov.subm.keys().collect();
            sids.sort();
            for sid in sids {
                let Some(_sb) = prov.subm.get(sid) else {
                    continue;
                };
                cn0 += 1;
                //let fnm = format!("/mnt/e/CHMBACK/pea-data/c01_pea/{sid}.bin");
                let Ok(buf) = std::fs::read(format!("{dnm}/{sid}.bin")) else {
                    continue;
                };
                //println!("      '{fnm}'");
                let (sub, _): (PeaSub, usize) =
                    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
                println!("      feed: {}", sub.feeders.len());
                cn1 += 1;

                /*
                let siv: Vec<_> = eg
                    .subs
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| v.sbid == *sid)
                    .map(|(i, _)| i)
                    .collect();
                if siv.len() != 1 {
                    continue;
                }
                let si = siv[0];
                */

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
                    let mut vf01 = 0f32;
                    let mut vf03 = f32::MAX;
                    let k1 = format!("{fid}");
                    let key = if let Some(f) = fd2fd.get(&k1) {
                        f.to_string()
                    } else {
                        "-".to_string()
                    };
                    if let Some(lp) = e0.lp24.get(&key) {
                        for v in lp.pos_rep.val.into_iter().flatten() {
                            vf01 = vf01.max(v.unwrap_or(0f32));
                            vf03 = vf03.min(v.unwrap_or(200f32));
                            if vf03 == f32::MAX {
                                vf03 = 0f32;
                            }
                        }
                    } else if let Some(lp) = e0.lp23.get(&key) {
                        for v in lp.pos_rep.val.into_iter().flatten() {
                            vf01 = vs01.max(v.unwrap_or(0f32));
                            vf03 = vf03.min(v.unwrap_or(200f32));
                            if vf03 == f32::MAX {
                                vf03 = 0f32;
                            }
                        }
                    };
                    vf03 = vf01 - vf03;

                    let mut vf02 = 0f32;
                    let mut vf04 = 0f32;
                    if let Some(lp) = e0.lp24.get(&key) {
                        for v in lp.neg_rep.val.into_iter().flatten() {
                            vf02 = vf02.max(v.unwrap_or(0f32));
                            vf04 = vf04.min(v.unwrap_or(200f32));
                            if vf04 == f32::MAX {
                                vf04 = 0f32;
                            }
                        }
                    } else if let Some(lp) = e0.lp23.get(&key) {
                        for v in lp.neg_rep.val.into_iter().flatten() {
                            vf02 = vs02.max(v.unwrap_or(0f32));
                            vf04 = vf04.min(v.unwrap_or(200f32));
                        }
                    };
                    vf04 = vf02 - vf04;

                    let mut tids: Vec<_> = fd.tranm.keys().collect();
                    tids.sort();
                    let mut tcn = 0;
                    let (mut mtc, mtz) = (0, 0);
                    for tid in tids {
                        let Some(trn) = fd.tranm.get(tid) else {
                            continue;
                        };
                        tcn += 1;
                        if trn.own == "P" {
                            txpc += 1;
                        } else {
                            txcc += 1;
                        }
                        mtc += trn.mets.len();
                        //let ev = trn.vopw;
                        //let kv = trn.tr_kva;
                        //let zi = trn.zons.clone();
                        let aojs = trn.aojs.clone();
                        let vt05 = trn.tr_kva.unwrap_or(10f32);
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
                                if met.main.is_empty() && met.kwh18 > 200f32 {
                                    vt10 += met.kwh15;
                                }
                            }
                            match met.mt_phs.clone().unwrap_or(String::new()).as_str() {
                                "A" => se_a += met.kwh15,
                                "B" => se_b += met.kwh15,
                                "C" => se_c += met.kwh15,
                                _ => {}
                            }
                        }
                        let vt11 = trn.mets.len() as f32;
                        let vt12 = 1f32;

                        let mut vt08 = 0f32;
                        let se_p = se_a + se_b + se_c;
                        if se_a < se_p && se_b < se_p && se_c < se_p {
                            let ab = (se_a - se_b).abs();
                            let bc = (se_b - se_c).abs();
                            let ca = (se_c - se_a).abs();
                            vt08 = (ab + bc + ca) * 0.5;
                        }
                        let vt09 = trf_kva_2_kw(vt02);

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

                        let tr_as = PeaAssVar {
                            arid: aid.to_string(),
                            pvid: pid.to_string(),
                            n1d: trn.n1d.to_string(),
                            sbid: sid.to_string(),
                            fdid: fid.to_string(),
                            own: trn.own.to_string(),
                            peano: trn.tr_pea.clone().unwrap_or("".to_string()).to_string(),
                            aoj,
                            vp01,
                            vp02,
                            vs01,
                            vs02,
                            vs03,
                            vs04,
                            vs05,
                            vs06,
                            vs07,
                            vf01,
                            vf02,
                            vf03,
                            vf04,
                            vt01,
                            vt02,
                            vt03,
                            vt04,
                            vt05,
                            vt06,
                            vt07,
                            vt08,
                            vt09,
                            vt10,
                            vt11,
                            vt12,
                            ..Default::default()
                        };
                        tras_mx1.max(&tr_as);
                        //s_tr_sum.add(&tr_as);
                        s_tr_ass.push(tr_as);
                    } // end trans loop
                    mtc0 += mtc;
                    mtz0 += mtz;
                    println!("      fd: {fid} tr:{tcn} mtc:{mtc} key:{k1}-{key}");
                } // end feeder loop
                write_trn_ass_01(&s_tr_ass, &format!("{dnm}/{sid}-raw.txt"))?;
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&s_tr_ass, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-raw.bin"), bin).unwrap();
            } // end sub loop
        } // end provi loop
    } // end area
    println!("c0:{cn0} c1:{cn1} {mtc0}-z:{mtz0} tp:{txpc} tc:{txcc}");

    ////////////////////////////////////////////////
    ////////////////////////////////////////////////
    let mut tras_mx2 = PeaAssVar::default();
    for id in &aids {
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
                //let fnm = format!("/mnt/e/CHMBACK/pea-data/c01_pea/{sid}-raw.bin");
                let Ok(buf) = std::fs::read(format!("{dnm}/{sid}-raw.bin")) else {
                    continue;
                };
                let (mut v_tras, _): (Vec<PeaAssVar>, usize) =
                    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
                println!("   {sid} - {}", v_tras.len());
                // normalize data
                for tras in &mut v_tras {
                    tras.nor(&tras_mx1);
                }
                //// save normal bin
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-nor.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras, &format!("{dnm}/{sid}-nor.txt"))?;

                ////////////////////////////////////////////////
                ////////////////////////////////////////////////
                // calculate EV
                let mut v_tras_ev = v_tras.clone();
                for (tras, tras0) in v_tras_ev.iter_mut().zip(v_tras.iter_mut()) {
                    tras.weigh(&evwe);
                    tras.sum();
                    tras0.vc01 = tras.sum;
                }
                //// save ev bin
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras_ev, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-ev.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras_ev, &format!("{dnm}/{sid}-ev.txt"))?;

                ////////////////////////////////////////////////
                ////////////////////////////////////////////////
                // calculate solar
                let mut v_tras_so = v_tras.clone();
                for (tras, tras0) in v_tras_so.iter_mut().zip(v_tras.iter_mut()) {
                    //for tras in &mut v_tras_so {
                    tras.weigh(&sowe);
                    tras.sum();
                    tras0.vc07 = tras.sum;
                }
                //// save ev bin
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras_ev, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-so.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras_ev, &format!("{dnm}/{sid}-so.txt"))?;

                for tr in v_tras.iter_mut() {
                    tr.vc02 = tr.vt09 / z2o(tr.vt05);
                    tr.vc03 = if tr.vc02 > 0.5f32 { 1f32 } else { 0f32 };
                    tr.vc04 = tr.vt03;
                    tr.vc05 = tr.vt04;
                    tr.vc06 = tr.vs01 / z2o(tr.vs07);
                    tr.vc08 = tr.vs03;
                    tr.vc09 = tr.vs04;
                    tr.vc10 = tr.vt02;
                    tr.vc11 = tr.vt10;
                    tr.vc12 = tr.vt08;
                    let v = tr.vt08 / z2o(tr.vt04);
                    tr.vc03 = if v > 0.5f32 { 1f32 } else { 0f32 };
                    tras_mx2.max(tr);
                }
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-rw2.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras, &format!("{dnm}/{sid}-rw2.txt"))?;
            } // end sub loop
        } // end provi loop
    } // end area

    for id in &aids {
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
                let (mut v_tras, _): (Vec<PeaAssVar>, usize) =
                    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
                println!("   {sid} - {}", v_tras.len());
                // normalize data
                for tras in &mut v_tras {
                    tras.nor(&tras_mx2);
                }
                //// save normal bin
                let bin: Vec<u8> =
                    bincode::encode_to_vec(&v_tras, bincode::config::standard()).unwrap();
                std::fs::write(format!("{dnm}/{sid}-no2.bin"), bin).unwrap();
                write_trn_ass_01(&v_tras, &format!("{dnm}/{sid}-no2.txt"))?;
            }
        }
    }

    Ok(())
}

fn write_trn_ass_01(tr_as: &Vec<PeaAssVar>, fnm: &str) -> Result<(), Box<dyn Error>> {
    let mut x = String::new();
    use std::fmt::Write;
    for t in tr_as {
        write!(x, "{}", t.sbid)?;
        write!(x, "\t{}", t.fdid)?;
        write!(x, "\t{}", t.aoj)?;
        write!(x, "\t{}", t.own)?;
        write!(x, "\t{}", t.peano)?;
        write!(x, "\t{}", t.vp01)?;
        write!(x, "\t{}", t.vp02)?;
        write!(x, "\t{}", t.vs01)?;
        write!(x, "\t{}", t.vs02)?;
        write!(x, "\t{}", t.vs03)?;
        write!(x, "\t{}", t.vs04)?;
        write!(x, "\t{}", t.vs05)?;
        write!(x, "\t{}", t.vs06)?;
        write!(x, "\t{}", t.vs07)?;
        write!(x, "\t{}", t.vf01)?;
        write!(x, "\t{}", t.vf02)?;
        write!(x, "\t{}", t.vf03)?;
        write!(x, "\t{}", t.vf04)?;
        write!(x, "\t{}", t.vt01)?;
        write!(x, "\t{}", t.vt02)?;
        write!(x, "\t{}", t.vt03)?;
        write!(x, "\t{}", t.vt04)?;
        write!(x, "\t{}", t.vt05)?;
        write!(x, "\t{}", t.vt06)?;
        write!(x, "\t{}", t.vt07)?;
        write!(x, "\t{}", t.vt08)?;
        write!(x, "\t{}", t.vt09)?;
        write!(x, "\t{}", t.vt10)?;
        write!(x, "\t{}", t.vt11)?;
        write!(x, "\t{}", t.vt12)?;
        write!(x, "\t{}", t.vc01)?;
        write!(x, "\t{}", t.vc02)?;
        write!(x, "\t{}", t.vc03)?;
        write!(x, "\t{}", t.vc04)?;
        write!(x, "\t{}", t.vc05)?;
        write!(x, "\t{}", t.vc06)?;
        write!(x, "\t{}", t.vc07)?;
        write!(x, "\t{}", t.vc08)?;
        write!(x, "\t{}", t.vc09)?;
        write!(x, "\t{}", t.vc10)?;
        write!(x, "\t{}", t.vc11)?;
        write!(x, "\t{}", t.vc12)?;
        write!(x, "\t{}", t.sum)?;
        writeln!(x)?;
    }
    println!("        ===== write to {fnm}");
    std::fs::write(fnm, x)?;
    Ok(())
}

fn write_trn_ass_02(tr_as: &Vec<PeaAssVar>, fnm: &str) -> Result<(), Box<dyn Error>> {
    let mut x = String::new();
    use std::fmt::Write;
    for t in tr_as {
        write!(x, "{}", t.sbid)?;
        write!(x, "\t{}", t.vp01)?;
        write!(x, "\t{}", t.vp02)?;
        write!(x, "\t{}", t.vs01)?;
        write!(x, "\t{}", t.vs02)?;
        write!(x, "\t{}", t.vs03)?;
        write!(x, "\t{}", t.vs04)?;
        write!(x, "\t{}", t.vs05)?;
        write!(x, "\t{}", t.vs06)?;
        write!(x, "\t{}", t.vs07)?;
        write!(x, "\t{}", t.vf01)?;
        write!(x, "\t{}", t.vf02)?;
        write!(x, "\t{}", t.vf03)?;
        write!(x, "\t{}", t.vf04)?;
        write!(x, "\t{}", t.vt01)?;
        write!(x, "\t{}", t.vt02)?;
        write!(x, "\t{}", t.vt03)?;
        write!(x, "\t{}", t.vt04)?;
        write!(x, "\t{}", t.vt05)?;
        write!(x, "\t{}", t.vt06)?;
        write!(x, "\t{}", t.vt07)?;
        write!(x, "\t{}", t.vt08)?;
        write!(x, "\t{}", t.vt09)?;
        write!(x, "\t{}", t.vt10)?;
        write!(x, "\t{}", t.vt11)?;
        write!(x, "\t{}", t.vt12)?;
        write!(x, "\t{}", t.vc01)?;
        write!(x, "\t{}", t.vc02)?;
        write!(x, "\t{}", t.vc03)?;
        write!(x, "\t{}", t.vc04)?;
        write!(x, "\t{}", t.vc05)?;
        write!(x, "\t{}", t.vc06)?;
        write!(x, "\t{}", t.vc07)?;
        write!(x, "\t{}", t.vc08)?;
        write!(x, "\t{}", t.vc09)?;
        write!(x, "\t{}", t.vc10)?;
        write!(x, "\t{}", t.vc11)?;
        write!(x, "\t{}", t.vc12)?;
        write!(x, "\t{}", t.sum)?;
        writeln!(x)?;
    }
    println!("        ===== write to {fnm}");
    std::fs::write(fnm, x)?;
    Ok(())
}
use std::collections::HashSet;

/// ประมวลผลรวมเพื่อเกณฑ์การคัดเลือก
/// summery transformaters to substation
pub fn c01_chk_03() -> Result<(), Box<dyn Error>> {
    let dnm = "/mnt/e/CHMBACK/pea-data/c01_pea";
    let buf = std::fs::read(format!("{dnm}/000_pea.bin")).unwrap();
    let (pea, _): (Pea, usize) =
        bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
    let mut aids: Vec<_> = pea.aream.keys().collect();
    aids.sort();
    let mut pvcn = 0;
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
            //println!("  pv:{pid}");
            pvcn += 1;
            let mut sids: Vec<_> = prov.subm.keys().collect();
            sids.sort();
            for sid in sids {
                let Some(_sb) = prov.subm.get(sid) else {
                    continue;
                };
                let Ok(buf) = std::fs::read(format!("{dnm}/{sid}-rw2.bin")) else {
                    continue;
                };
                let (v_tras, _): (Vec<PeaAssVar>, usize) =
                    bincode::decode_from_slice(&buf[..], bincode::config::standard()).unwrap();
                if v_tras.is_empty() {
                    println!("    {sid} - NO data ");
                    continue;
                }
                let tras = &v_tras[0];
                let mut sbas = PeaAssVar {
                    arid: aid.to_string(),
                    pvid: pid.to_string(),
                    sbid: tras.sbid.to_string(),
                    ..Default::default()
                };
                for tras in &v_tras {
                    sbas.add(tras);
                }
                sbas.vp01 = tras.vp01;
                sbas.vp02 = tras.vp02;
                sbas.vs01 = tras.vs01;
                sbas.vs02 = tras.vs02;
                sbas.vs03 = tras.vs03;
                sbas.vs04 = tras.vs04;
                sbas.vs05 = tras.vs05;
                sbas.vs06 = tras.vs06;
                sbas.vs07 = tras.vs07;
                sbas_mx.max(&sbas);
                v_sbas.push(sbas);
                //println!("   {sid} - {}", v_tras.len());
            } // end sub loop
        } // end provi loop
    } // end area
    println!("pvcn: {pvcn}");

    // save ev bin
    let bin: Vec<u8> = bincode::encode_to_vec(&v_sbas, bincode::config::standard()).unwrap();
    std::fs::write(format!("{dnm}/000-sbrw.bin"), bin).unwrap();
    write_trn_ass_02(&v_sbas, &format!("{dnm}/000-sbrw.txt"))?;

    for sub in v_sbas.iter_mut() {
        sub.nor(&sbas_mx);
    }
    let bin: Vec<u8> = bincode::encode_to_vec(&v_sbas, bincode::config::standard()).unwrap();
    std::fs::write(format!("{dnm}/000-sbno.bin"), bin).unwrap();
    write_trn_ass_02(&v_sbas, &format!("{dnm}/000-sbno.txt"))?;

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
