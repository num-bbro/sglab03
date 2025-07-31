use crate::p02::DrawLoadProf;
use crate::p02::LoadProf;
use crate::p02::SubLoadProf;
use crate::p03::SubLoadProfRepr;
use std::collections::HashMap;
use std::error::Error;

use bincode::Decode;
use bincode::Encode;
use strum_macros::EnumIter;

use iterstats::Iterstats;
use num::complex::Complex;
use regex::Regex;
use std::f64::consts::PI;
const FK1_RT_SOLA: f32 = 0.15;
const BIGGER_RATIO: f32 = 1.5;

pub fn p08_draw_01(yr: &str, _sb: &str) -> Result<(), Box<dyn Error>> {
    let fnm = format!("/mnt/e/CHMBACK/pea-data/data2/p03_lp_repr_{yr}.bin");
    println!("{yr} {fnm}");
    let bytes = std::fs::read(fnm).unwrap();
    let (sblpr, _): (HashMap<String, SubLoadProfRepr>, usize) =
        bincode::decode_from_slice(&bytes[..], bincode::config::standard()).unwrap();
    println!("sblpr:{}", sblpr.len());
    //let re1 = Regex::new(r"^[A-Z]{3}[0-1][0-9][VW][BR]01$").unwrap();
    let re1 = Regex::new(r"^[A-Z]{3}-[0-1][0-9A-Z][VW][BR]01$").unwrap();
    let re2 = Regex::new(r"^[A-Z]{3}$").unwrap();
    for (i, (sb, slp)) in sblpr.iter().enumerate() {
        if !re1.is_match(sb) && !re2.is_match(sb) {
            continue;
        }
        if let Some(lpr) = &slp.pos_rep.val {
            use std::fmt::Write;
            let mut ss = String::new();
            print!("{i}. {sb}");
            //let mut lpr = lpr.clone();
            let mut lpr = *lpr;
            for i in 0..lpr.len() {
                if lpr[i].is_none() {
                    if i == 0 {
                        print!("_");
                        lpr[i] = lpr[1];
                    } else if i == lpr.len() - 1 {
                        print!(".");
                        lpr[i] = lpr[lpr.len() - 2];
                    } else if let (Some(n1), Some(n2)) = (lpr[i - 1], lpr[i + 1]) {
                        print!("=");
                        lpr[i] = Some((n1 + n2) * 0.5f32);
                    }
                }
            }
            let nc: Vec<_> = lpr.iter().flatten().map(|v| *v).collect();
            //let nc: Vec<_> = lpr.iter().flatten().copied();
            let mx = lpr.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            println!(" x:{mx:?}");
            if nc.len() == lpr.len() {
                print!("    F:");
                #[allow(non_camel_case_types)]
                let nn = nc.len();
                let (mut en, mut en1, mut en2) = (0f32, 0f32, 0f32);
                const DFT_K_MAX: usize = 4;
                let mut ffk = [Complex::new(0.0, 0.0); DFT_K_MAX];
                for (k, ffk) in ffk.iter_mut().enumerate() {
                    print!(" Xk{k}:");
                    let mut xxk = Complex::new(0f64, 0f64);
                    for (n, v) in nc.iter().enumerate() {
                        en += *v;
                        en1 += if *v < 0f32 { -*v } else { 0f32 };
                        en2 += if *v > 0f32 { *v } else { 0f32 };
                        let xn = Complex::new(*v as f64, 0f64);
                        let ex = Complex::new(0f64, -2f64 * PI * k as f64 * n as f64 / nn as f64);
                        let ex = xn * ex.exp();
                        xxk += ex;
                        //print!("{ex} ");
                    }
                    *ffk = xxk / nn as f64;
                    //println!();
                    print!("{ffk:.2} ");
                }
                println!();
                let min = nc.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
                let max = nc.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
                let dif = max - min;
                let fk1 = ffk[1].re.abs() / dif as f64;
                let pks = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i > 24 && *i < 72)
                    .map(|(_, v)| -v)
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();
                let ngm = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i < 24)
                    .map(|(_, v)| v)
                    .mean();
                let eno: f32 = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i >= 72)
                    .map(|(_, v)| v - ngm)
                    .sum();
                let enb = ngm * 6.0;
                let eno = eno * 0.25;
                let bor = eno / enb;
                let pt1: f32 = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i < 24)
                    .map(|(_, v)| v)
                    .sum();
                let pt2: f32 = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i >= 24 && *i < 48)
                    .map(|(_, v)| v)
                    .sum();
                let pt3: f32 = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i >= 48 && *i < 72)
                    .map(|(_, v)| v)
                    .sum();
                let pt4: f32 = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i >= 72)
                    .map(|(_, v)| v)
                    .sum();
                //
                //let mut rc = [Complex::new(0f64, 0f64); 96];
                //print!("    f:");
                //println!();

                let mut gp = "GZ";
                let po = ffk[1].to_polar();
                if po.0.abs() + ffk[0].re.abs() < 0.1 {
                    gp = "G0";
                } else if en < 0.0 && en2 < 0.1 * en1 {
                    gp = "G1";
                }
                let ang = po.1 / PI * 180f64;
                //let rt1 = ffk[1].re / ffk[0].re;
                if en2 > 0.0 && fk1 >= FK1_RT_SOLA.into() && ang > -70.0 && ang < 30.0 {
                    //if en2 > 0.0 && ffk[1].re > 0.1 && ang > -70.0 && ang < 30.0 {
                    //if rt1 > 0.2 && en2 > 0.0 && ffk[1].re > 0.1 && ang > -70.0 && ang < 30.0 {
                    gp = "G2";
                }
                if gp == "G2" && enb > 0.0 && bor > 0.3 {
                    gp = "G3";
                }
                if gp == "GZ" && enb > 0.0 && bor > 0.5 && (pt3 + pt4) / (pt1 + pt2) > BIGGER_RATIO
                {
                    gp = "G4";
                }
                let pt0 = pt1 > 0.0 && pt2 > 0.0 && pt3 > 0.0 && pt4 > 0.0;
                if gp == "GZ" && pt0 && (pt1 + pt2) / (pt3 + pt4) > BIGGER_RATIO {
                    gp = "G5";
                }
                if gp == "GZ" && pt0 && (pt3 + pt4) / (pt1 + pt2) > BIGGER_RATIO {
                    gp = "G6";
                }
                if gp == "GZ" && pt0 && (pt2 + pt3) / (pt1 + pt4) > BIGGER_RATIO {
                    gp = "G7";
                }
                if gp == "GZ" && pt0 && (pt1 + pt4) / (pt2 + pt3) > BIGGER_RATIO {
                    gp = "G8";
                }
                if gp == "GZ" && pt1 < 0.0 && pt2 > 0.0 && pt3 > 0.0 && pt4 < 0.0 {
                    gp = "G9";
                }
                if gp == "GZ" && pt1 < 0.0 && pt2 < 0.0 && pt3 < 0.0 && pt4 > 0.0 {
                    gp = "GA";
                }
                if gp == "G1" && fk1 >= FK1_RT_SOLA.into() && ang < -150.0 && ang > -210.0 {
                    gp = "GB";
                }
                let mut rc = [0f64; 96];
                if gp == "G2" || gp == "G3" {
                    for (n, v) in rc.iter_mut().enumerate() {
                        if n < 24 || n > 72 {
                            continue;
                        }
                        let mut vn = Complex::new(0f64, 0f64);
                        for (k, ffk) in ffk.iter().enumerate() {
                            //if k != 1 {
                            if k < 1 {
                                continue;
                            }
                            let ex =
                                Complex::new(0f64, 2f64 * PI * k as f64 * n as f64 / nn as f64);
                            let mut ffk0 = *ffk;
                            if k == 1 {
                                let fk = ffk.to_polar().0;
                                let mu = (pks + ngm) * 0.9;
                                ffk0 = ffk0 / fk * mu as f64;
                            }
                            let ex = ffk0 * ex.exp();
                            vn += ex;
                        }
                        if vn.re < 0.0 {
                            *v = vn.re;
                        }
                        //let re = vn.re;
                        //print!("{re:.2}[{:.2}] ", ffk[0].re);
                    }
                }

                writeln!(ss, "all energy: {en}")?;
                writeln!(ss, "neg energy: {en1}")?;
                writeln!(ss, "pos energy: {en2}")?;
                writeln!(ss, "max power: {max}")?;
                writeln!(ss, "min power: {min}")?;
                writeln!(ss, "dif power: {dif}")?;

                writeln!(ss, "Fk0 Re = {:.2}", ffk[0].re)?;
                writeln!(ss, "Fk1 rt = {:.2}", fk1)?;
                writeln!(ss, "Fk0 Im = {:.2}", ffk[0].im)?;
                writeln!(ss, "Fk1 Re = {:.2}", ffk[1].re)?;
                writeln!(ss, "Fk1 Im = {:.2}", ffk[1].im)?;
                writeln!(ss, "Fk1 Mag = {:.2}", po.0)?;
                writeln!(ss, "Fk1 Ang = {:.2}", po.1 / PI * 180f64)?;
                //let ngm = nc.iter().mean();
                writeln!(ss, "Night avg: {:.2} mw", ngm)?;
                writeln!(ss, "Base eng: {:.2} mwh", enb)?;
                writeln!(ss, "Over eng: {:.2} mwh", eno)?;
                writeln!(ss, "Over/Base: {:.2} %", bor)?;
                writeln!(ss, "Peak sola: {:.2} mw", pks)?;

                let dnm = format!("/mnt/e/CHMBACK/pea-data/p08_drw_po/{yr}/{gp}-{sb}");
                std::fs::create_dir_all(&dnm)?;

                let fnm = format!("{dnm}/para.txt");
                std::fs::write(fnm, ss)?;

                ///////////////////////////////////////
                // NEG1
                let mut val = [Some(0f32); 96];
                for (i, v) in nc.iter().enumerate() {
                    val[i] = Some(*v);
                }
                let lb1 = format!("{sb}-FEQ0");
                let lb2 = "NEG".to_string();
                let fnm = format!("{dnm}/{lb1}.jpg");
                println!("WR -> {fnm}");
                let sz = (400, 300);
                let rf = Vec::<(String, f32)>::new();

                let drlp = LoadProf {
                    lb1,
                    lb2,
                    fnm,
                    val,
                    sz,
                    rf,
                };
                let Ok(_) = drlp.draw_prof() else {
                    continue;
                };

                if gp == "G2" || gp == "G3" {
                    ///////////////////////////////////////
                    // NEG2
                    let mut val = [Some(0f32); 96];
                    for (i, v) in rc.iter().enumerate() {
                        val[i] = Some(*v as f32);
                    }
                    let lb1 = format!("{sb}-FEQ1");
                    let lb2 = "NEG".to_string();
                    //let val: Vec<_> = nc.iter().map(|a| Some(a)).collect();
                    let fnm = format!("{dnm}/{lb1}.jpg");
                    println!("WR -> {fnm}");
                    let sz = (400, 300);
                    let rf = Vec::<(String, f32)>::new();

                    let drlp = LoadProf {
                        lb1,
                        lb2,
                        fnm,
                        val,
                        sz,
                        rf,
                    };
                    let Ok(_) = drlp.draw_prof() else {
                        continue;
                    };
                }
            }
            //print!(" n:{}", nc.len());
        };
        /*
        if let Some(lpr) = &slp.neg_rep.val {
            use std::fmt::Write;
            let mut ss = String::new();
            print!("{i}. {sb}");
            //let mut lpr = lpr.clone();
            let mut lpr = *lpr;
            for i in 0..lpr.len() {
                if lpr[i].is_none() {
                    if i == 0 {
                        print!("_");
                        lpr[i] = lpr[1];
                    } else if i == lpr.len() - 1 {
                        print!(".");
                        lpr[i] = lpr[lpr.len() - 2];
                    } else if let (Some(n1), Some(n2)) = (lpr[i - 1], lpr[i + 1]) {
                        print!("=");
                        lpr[i] = Some((n1 + n2) * 0.5f32);
                    }
                }
            }
            let nc: Vec<_> = lpr.iter().flatten().map(|v| *v).collect();
            //let nc: Vec<_> = lpr.iter().flatten().copied();
            let mx = lpr.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            println!(" x:{mx:?}");
            if nc.len() == lpr.len() {
                print!("    F:");
                #[allow(non_camel_case_types)]
                let nn = nc.len();
                let (mut en, mut en1, mut en2) = (0f32, 0f32, 0f32);
                const DFT_K_MAX: usize = 4;
                let mut ffk = [Complex::new(0.0, 0.0); DFT_K_MAX];
                for (k, ffk) in ffk.iter_mut().enumerate() {
                    print!(" Xk{k}:");
                    let mut xxk = Complex::new(0f64, 0f64);
                    for (n, v) in nc.iter().enumerate() {
                        en += *v;
                        en1 += if *v < 0f32 { -*v } else { 0f32 };
                        en2 += if *v > 0f32 { *v } else { 0f32 };
                        let xn = Complex::new(*v as f64, 0f64);
                        let ex = Complex::new(0f64, -2f64 * PI * k as f64 * n as f64 / nn as f64);
                        let ex = xn * ex.exp();
                        xxk += ex;
                        //print!("{ex} ");
                    }
                    *ffk = xxk / nn as f64;
                    //println!();
                    print!("{ffk:.2} ");
                }
                println!();
                let pks = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i > 24 && *i < 72)
                    .map(|(_, v)| -v)
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();
                let ngm = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i < 24)
                    .map(|(_, v)| v)
                    .mean();
                let eno: f32 = nc
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i >= 72)
                    .map(|(_, v)| v - ngm)
                    .sum();
                let enb = ngm * 6.0;
                let eno = eno * 0.25;
                let bor = eno / enb;
                //
                //let mut rc = [Complex::new(0f64, 0f64); 96];
                //print!("    f:");
                //println!();
                let mut gp = "GZ";
                let po = ffk[1].to_polar();
                if po.0.abs() + ffk[0].re.abs() < 0.1 {
                    gp = "G0";
                } else if en < 0.0 && en2 < 0.1 * en1 {
                    gp = "G1";
                }
                let ang = po.1 / PI * 180f64;
                if en2 > 0.0 && ffk[1].re > 0.1 && ang > -70.0 && ang < 30.0 {
                    gp = "G2";
                }
                if gp == "G2" && enb > 0.0 && bor > 0.3 {
                    gp = "G3";
                }
                if gp == "GZ" && enb > 0.0 && bor > 0.5 {
                    gp = "G4";
                }
                let mut rc = [0f64; 96];
                if gp == "G2" || gp == "G3" {
                    for (n, v) in rc.iter_mut().enumerate() {
                        if n < 24 || n > 72 {
                            continue;
                        }
                        let mut vn = Complex::new(0f64, 0f64);
                        for (k, ffk) in ffk.iter().enumerate() {
                            //if k != 1 {
                            if k < 1 {
                                continue;
                            }
                            let ex =
                                Complex::new(0f64, 2f64 * PI * k as f64 * n as f64 / nn as f64);
                            let mut ffk0 = *ffk;
                            if k == 1 {
                                let fk = ffk.to_polar().0;
                                let mu = (pks + ngm) * 0.9;
                                ffk0 = ffk0 / fk * mu as f64;
                            }
                            let ex = ffk0 * ex.exp();
                            vn += ex;
                        }
                        if vn.re < 0.0 {
                            *v = vn.re;
                        }
                        //let re = vn.re;
                        //print!("{re:.2}[{:.2}] ", ffk[0].re);
                    }
                }

                writeln!(ss, "all energy: {en}")?;
                writeln!(ss, "neg energy: {en1}")?;
                writeln!(ss, "pos energy: {en2}")?;

                writeln!(ss, "Fk0 Re = {:.2}", ffk[0].re)?;
                writeln!(ss, "Fk0 Im = {:.2}", ffk[0].im)?;
                writeln!(ss, "Fk1 Re = {:.2}", ffk[1].re)?;
                writeln!(ss, "Fk1 Im = {:.2}", ffk[1].im)?;
                writeln!(ss, "Fk1 Mag = {:.2}", po.0)?;
                writeln!(ss, "Fk1 Ang = {:.2}", po.1 / PI * 180f64)?;
                //let ngm = nc.iter().mean();
                writeln!(ss, "Night avg: {:.2} mw", ngm)?;
                writeln!(ss, "Base eng: {:.2} mwh", enb)?;
                writeln!(ss, "Over eng: {:.2} mwh", eno)?;
                writeln!(ss, "Over/Base: {:.2} %", bor)?;
                writeln!(ss, "Peak sola: {:.2} mw", pks)?;

                let dnm = format!("/mnt/e/CHMBACK/pea-data/p08_drw_01/{yr}/{gp}-{sb}");
                std::fs::create_dir_all(&dnm)?;

                let fnm = format!("{dnm}/para.txt");
                std::fs::write(fnm, ss)?;
                ///////////////////////////////////////
                // NEG1
                let mut val = [Some(0f32); 96];
                for (i, v) in nc.iter().enumerate() {
                    val[i] = Some(*v);
                }
                let lb1 = format!("{sb}-FEQ0");
                let lb2 = "NEG".to_string();
                let fnm = format!("{dnm}/{lb1}.jpg");
                println!("WR -> {fnm}");
                let sz = (400, 300);
                let rf = Vec::<(String, f32)>::new();

                let drlp = LoadProf {
                    lb1,
                    lb2,
                    fnm,
                    val,
                    sz,
                    rf,
                };
                let Ok(_) = drlp.draw_prof() else {
                    continue;
                };

                ///////////////////////////////////////
                // NEG2
                let mut val = [Some(0f32); 96];
                for (i, v) in rc.iter().enumerate() {
                    val[i] = Some(*v as f32);
                }
                let lb1 = format!("{sb}-FEQ1");
                let lb2 = "NEG".to_string();
                //let val: Vec<_> = nc.iter().map(|a| Some(a)).collect();
                let fnm = format!("{dnm}/{lb1}.jpg");
                println!("WR -> {fnm}");
                let sz = (400, 300);
                let rf = Vec::<(String, f32)>::new();

                let drlp = LoadProf {
                    lb1,
                    lb2,
                    fnm,
                    val,
                    sz,
                    rf,
                };
                let Ok(_) = drlp.draw_prof() else {
                    continue;
                };
            }
            //print!(" n:{}", nc.len());
        };
        */
    }
    Ok(())
}

pub fn p08_draw_slp(yr: &str, sb: &str) -> Result<(), Box<dyn Error>> {
    let fnm = format!("/mnt/e/CHMBACK/pea-data/data2/p02_sub_lp_{yr}.bin");
    let bytes = std::fs::read(fnm).unwrap();
    let (subh, _): (HashMap<String, SubLoadProf>, usize) =
        bincode::decode_from_slice(&bytes[..], bincode::config::standard()).unwrap();
    println!("sub {}", subh.len());
    let mut keys: Vec<_> = subh.keys().collect();
    keys.sort();
    {
        let sb = sb.to_string();
        let sub = subh.get(&sb).unwrap();
        println!("s: {}", sub.sub);
        for f in &sub.feeds {
            let f = f.to_string();
            let flp = sub.fdldp.get(&f).unwrap();
            println!("  f:{f} - {}", flp.days.len());
            for (di, dlp) in flp.days.iter().enumerate() {
                if let Some(dlp) = dlp {
                    let lb1 = format!("{sb}-{f}");
                    let mon = dlp.mon;
                    let mdt = dlp.mdt;
                    let lb2 = format!("{yr}-{mon:02}-{mdt:02}");
                    let val = dlp.val;
                    let dnm = format!(
                        "/mnt/e/CHMBACK/pea-data/p08drw/{yr}/{sb}/{di:03}_{mon:02}-{mdt:02}"
                    );
                    std::fs::create_dir_all(&dnm)?;
                    let fnm = format!("{dnm}/{lb1}.jpg");
                    println!("{lb1} {lb2} -> {fnm}");
                    let sz = (400, 300);
                    let rf = Vec::<(String, f32)>::new();

                    let drlp = LoadProf {
                        lb1,
                        lb2,
                        fnm,
                        val,
                        sz,
                        rf,
                    };
                    let Ok(_) = drlp.draw_prof() else {
                        continue;
                    };
                }
            }
        }

        if let Some(ldpf) = &sub.ldpf {
            for (di, dlp) in ldpf.days.iter().enumerate() {
                if let Some(dlp) = dlp {
                    let lb1 = sb.to_string();
                    let mon = dlp.mon;
                    let mdt = dlp.mdt;
                    let lb2 = format!("{yr}-{mon:02}-{mdt:02}");
                    let val = dlp.val;
                    let dnm = format!(
                        "/mnt/e/CHMBACK/pea-data/p08drw/{yr}/{sb}/{di:03}_{mon:02}-{mdt:02}"
                    );
                    std::fs::create_dir_all(&dnm)?;
                    let fnm = format!("{dnm}/{lb1}.jpg");
                    println!("{lb1} {lb2} -> {fnm}");
                    let sz = (400, 300);
                    let rf = Vec::<(String, f32)>::new();

                    let drlp = LoadProf {
                        lb1,
                        lb2,
                        fnm,
                        val,
                        sz,
                        rf,
                    };
                    let Ok(_) = drlp.draw_prof() else {
                        continue;
                    };
                }
            }
        } else {
            println!("  NO SUB SUM");
        }
    }

    Ok(())
}

use crate::p02::DAY_VAL_PNTS;
#[derive(Encode, Decode, Debug, Clone, Default, EnumIter, Eq, Hash, PartialEq)]
pub enum ProfType {
    #[default]
    NoPower, // G0
    Injection,    // G1
    SolarPower,   // G2
    SolarNight,   // G3
    Night,        // G4
    FirstHalf,    // G5
    LastHalf,     // G6
    MiddHalf,     // G7
    BothEndHalf,  // G8
    NegPosPosNeg, // G9
    NegNegPosPos, // GA
    InjectSolar,  // GB
    Unknown,      // GZ
}

#[derive(Encode, Decode, Debug, Clone, Default)]
pub struct ProfInfo {
    feed: String,
    val: Vec<f32>,
    lp_type: ProfType,
    all_en: f32,
    pos_en: f32,
    neg_en: f32,
    max_pw: f32,
    min_pw: f32,
    dif_pw: f32,
    fk0_re: f32,
    fk1_re: f32,
    fk1_mg: f32,
    fk1_an: f32,
    ngt_av: f32,
    bas_en: f32,
    ovr_en: f32,
    ovr_rt: f32,
    sol_pk: Option<f32>,
    sol_en: Option<f32>,
    solar: Option<Vec<f32>>,
}

pub fn p08_class_val(lpv: &[Option<f32>; DAY_VAL_PNTS]) -> Result<ProfInfo, Box<dyn Error>> {
    let mut lpr = *lpv;
    for i in 0..lpr.len() {
        if lpr[i].is_none() {
            if i == 0 {
                lpr[i] = lpr[1];
            } else if i == lpr.len() - 1 {
                lpr[i] = lpr[lpr.len() - 2];
            } else if let (Some(n1), Some(n2)) = (lpr[i - 1], lpr[i + 1]) {
                lpr[i] = Some((n1 + n2) * 0.5f32);
            }
        }
    }
    let nc: Vec<_> = lpr.iter().flatten().map(|v| *v).collect();
    let mx = lpr.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    //println!(" x:{mx:?}");
    if nc.len() == lpr.len() {
        //print!("    F:");
        #[allow(non_camel_case_types)]
        let nn = nc.len();
        let (mut en, mut en1, mut en2) = (0f32, 0f32, 0f32);
        const DFT_K_MAX: usize = 4;
        let mut ffk = [Complex::new(0.0, 0.0); DFT_K_MAX];
        for (k, ffk) in ffk.iter_mut().enumerate() {
            //print!(" Xk{k}:");
            let mut xxk = Complex::new(0f64, 0f64);
            for (n, v) in nc.iter().enumerate() {
                en += *v;
                en1 += if *v < 0f32 { -*v } else { 0f32 };
                en2 += if *v > 0f32 { *v } else { 0f32 };
                let xn = Complex::new(*v as f64, 0f64);
                let ex = Complex::new(0f64, -2f64 * PI * k as f64 * n as f64 / nn as f64);
                let ex = xn * ex.exp();
                xxk += ex;
                //print!("{ex} ");
            }
            *ffk = xxk / nn as f64;
            //println!();
            //print!("{ffk:.2} ");
        }
        //println!();
        let min = nc.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let max = nc.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        let dif = max - min;
        let fk1 = ffk[1].re.abs() / dif as f64;
        let pks = nc
            .iter()
            .enumerate()
            .filter(|(i, _)| *i > 24 && *i < 72)
            .map(|(_, v)| -v)
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let ngm = nc
            .iter()
            .enumerate()
            .filter(|(i, _)| *i < 24)
            .map(|(_, v)| v)
            .mean();
        let eno: f32 = nc
            .iter()
            .enumerate()
            .filter(|(i, _)| *i >= 72)
            .map(|(_, v)| v - ngm)
            .sum();
        let enb = ngm * 6.0;
        let eno = eno * 0.25;
        let bor = eno / enb;
        let pt1: f32 = nc
            .iter()
            .enumerate()
            .filter(|(i, _)| *i < 24)
            .map(|(_, v)| v)
            .sum();
        let pt2: f32 = nc
            .iter()
            .enumerate()
            .filter(|(i, _)| *i >= 24 && *i < 48)
            .map(|(_, v)| v)
            .sum();
        let pt3: f32 = nc
            .iter()
            .enumerate()
            .filter(|(i, _)| *i >= 48 && *i < 72)
            .map(|(_, v)| v)
            .sum();
        let pt4: f32 = nc
            .iter()
            .enumerate()
            .filter(|(i, _)| *i >= 72)
            .map(|(_, v)| v)
            .sum();
        //
        let mut lptp = ProfType::Unknown;
        let po = ffk[1].to_polar();
        if po.0.abs() + ffk[0].re.abs() < 0.1 {
            lptp = ProfType::NoPower;
        } else if en < 0.0 && en2 < 0.1 * en1 {
            lptp = ProfType::Injection;
        }
        let ang = po.1 / PI * 180f64;
        if en2 > 0.0 && fk1 >= FK1_RT_SOLA.into() && ang > -70.0 && ang < 30.0 {
            lptp = ProfType::SolarPower;
        }
        if lptp == ProfType::SolarPower && enb > 0.0 && bor > 0.3 {
            lptp = ProfType::SolarNight;
        }
        if lptp == ProfType::Unknown
            && enb > 0.0
            && bor > 0.5
            && (pt3 + pt4) / (pt1 + pt2) > BIGGER_RATIO
        {
            lptp = ProfType::Night;
        }
        let pt0 = pt1 > 0.0 && pt2 > 0.0 && pt3 > 0.0 && pt4 > 0.0;
        if lptp == ProfType::Unknown && pt0 && (pt1 + pt2) / (pt3 + pt4) > BIGGER_RATIO {
            lptp = ProfType::FirstHalf;
        }
        if lptp == ProfType::Unknown && pt0 && (pt3 + pt4) / (pt1 + pt2) > BIGGER_RATIO {
            lptp = ProfType::LastHalf;
        }
        if lptp == ProfType::Unknown && pt0 && (pt2 + pt3) / (pt1 + pt4) > BIGGER_RATIO {
            lptp = ProfType::MiddHalf;
        }
        if lptp == ProfType::Unknown && pt0 && (pt1 + pt4) / (pt2 + pt3) > BIGGER_RATIO {
            lptp = ProfType::BothEndHalf;
        }
        if lptp == ProfType::Unknown && pt1 < 0.0 && pt2 > 0.0 && pt3 > 0.0 && pt4 < 0.0 {
            lptp = ProfType::NegPosPosNeg;
        }
        if lptp == ProfType::Unknown && pt1 < 0.0 && pt2 < 0.0 && pt3 < 0.0 && pt4 > 0.0 {
            lptp = ProfType::NegNegPosPos;
        }
        if lptp == ProfType::Injection && fk1 >= FK1_RT_SOLA.into() && ang < -150.0 && ang > -210.0
        {
            lptp = ProfType::InjectSolar;
        }
        let mut rc = [0f64; 96];
        //if [ProfType::SolarPower, ProfType::SolarNight].contains(lptp) {
        if lptp == ProfType::SolarPower || lptp == ProfType::SolarNight {
            for (n, v) in rc.iter_mut().enumerate() {
                if n < 24 || n > 72 {
                    continue;
                }
                let mut vn = Complex::new(0f64, 0f64);
                for (k, ffk) in ffk.iter().enumerate() {
                    //if k != 1 {
                    if k < 1 {
                        continue;
                    }
                    let ex = Complex::new(0f64, 2f64 * PI * k as f64 * n as f64 / nn as f64);
                    let mut ffk0 = *ffk;
                    if k == 1 {
                        let fk = ffk.to_polar().0;
                        let mu = (pks + ngm) * 0.9;
                        ffk0 = ffk0 / fk * mu as f64;
                    }
                    let ex = ffk0 * ex.exp();
                    vn += ex;
                }
                if vn.re < 0.0 {
                    *v = vn.re;
                }
                //let re = vn.re;
                //print!("{re:.2}[{:.2}] ", ffk[0].re);
            }
        }
        let (mut sol_pk, mut sol_en, mut solar) = (None, None, None);
        if lptp == ProfType::SolarPower || lptp == ProfType::SolarNight {
            let so: Vec<_> = rc.iter().map(|v| *v as f32).collect();
            let pk = so.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
            let pk = -*pk;
            sol_pk = Some(pk);
            let en: f32 = so.iter().sum();
            sol_en = Some(en);
            solar = Some(so);
        }
        let lpinf = ProfInfo {
            lp_type: lptp,
            val: nc.clone(),
            all_en: en,
            pos_en: en1,
            neg_en: en2,
            max_pw: *max,
            min_pw: *min,
            dif_pw: dif,
            fk0_re: ffk[0].re as f32,
            fk1_re: ffk[1].re as f32,
            fk1_mg: po.0 as f32,
            fk1_an: po.1 as f32,
            ngt_av: ngm,
            bas_en: enb,
            ovr_en: eno,
            ovr_rt: bor,
            sol_pk,
            sol_en,
            solar,
            ..Default::default()
        };
        return Ok(lpinf);
    }
    Err("Error".into())
}

pub fn p08_calc_lp1(yr: &str) -> Result<(), Box<dyn Error>> {
    let fnm = format!("/mnt/e/CHMBACK/pea-data/data2/p02_sub_lp_{yr}.bin");
    let bytes = std::fs::read(fnm).unwrap();
    let (subh, _): (HashMap<String, SubLoadProf>, usize) =
        bincode::decode_from_slice(&bytes[..], bincode::config::standard()).unwrap();
    //let mut mon_fdldp = HashMap::<usize,F
    let (mut n1w, mut n2w) = (0, 0);
    let (mut n1e, mut n2e) = (0.0, 0.0);
    for (s, sub) in subh {
        //println!("{s}");
        if let Some(_ldpf) = &sub.ldpf {
            //
        }
        for ldpf in sub.fdldp.values() {
            let f = &ldpf.feed;
            let mut lpv345 = Vec::<ProfInfo>::new();
            let mut lpv678 = Vec::<ProfInfo>::new();
            for (di, dlp) in ldpf.days.iter().enumerate() {
                if let Some(dlp) = dlp {
                    if let Ok(lpf) = p08_class_val(&dlp.val) {
                        if lpf.lp_type == ProfType::SolarPower
                            || lpf.lp_type == ProfType::SolarNight
                        {
                            /////////////////////////////////////
                            let mut val = [None; 96];
                            let lb1 = format!("{s}-{f}");
                            let mon = dlp.mon;
                            let mdt = dlp.mdt;
                            let lb2 = format!("{yr}-{mon:02}-{mdt:02}");
                            for (i, v) in lpf.val.iter().enumerate() {
                                val[i] = Some(*v);
                            }
                            let dnm = format!( "/mnt/e/CHMBACK/pea-data/p08sol/{yr}/{s}-{f}/{di:03}_{mon:02}-{mdt:02}");
                            std::fs::create_dir_all(&dnm)?;
                            let fnm = format!("{dnm}/prof.jpg");
                            println!("{lb1} {lb2} -> {fnm}");
                            let sz = (400, 300);
                            let rf = Vec::<(String, f32)>::new();
                            let drlp = LoadProf {
                                lb1,
                                lb2,
                                fnm,
                                val,
                                sz,
                                rf,
                            };
                            let Ok(_) = drlp.draw_prof() else {
                                continue;
                            };

                            /////////////////////////////////////
                            let mut val = [None; 96];
                            let lb1 = format!("{s}-{f} Solar");
                            let mon = dlp.mon;
                            let mdt = dlp.mdt;
                            let lb2 = format!("{yr}-{mon:02}-{mdt:02}");
                            let Some(solar) = &lpf.solar else {
                                continue;
                            };
                            for (i, v) in solar.iter().enumerate() {
                                val[i] = Some(*v);
                            }
                            let dnm = format!( "/mnt/e/CHMBACK/pea-data/p08sol/{yr}/{s}-{f}/{di:03}_{mon:02}-{mdt:02}");
                            std::fs::create_dir_all(&dnm)?;
                            let fnm = format!("{dnm}/solar.jpg");
                            println!("{lb1} {lb2} -> {fnm}");
                            let sz = (400, 300);
                            let rf = Vec::<(String, f32)>::new();
                            let drlp = LoadProf {
                                lb1,
                                lb2,
                                fnm,
                                val,
                                sz,
                                rf,
                            };
                            let Ok(_) = drlp.draw_prof() else {
                                continue;
                            };

                            if [4, 5].contains(&dlp.mon) {
                                lpv345.push(lpf);
                            } else if [6, 7].contains(&dlp.mon) {
                                lpv678.push(lpf);
                            }
                        }
                        //println!("    {i}. m:{} d:{} tp:{:?}", dlp.mon, dlp.mdt, lpf.lp_type);
                    } else {
                        //println!("    {i}. Error m:{} d:{}", dlp.mon, dlp.mdt);
                    }
                } else {
                    //println!("    {i}. ============");
                }
            }
            if lpv345.len() > 5 && lpv678.len() > 5 {
                let n1 = lpv345.len() as f32;
                let n2 = lpv678.len() as f32;
                let e1: f32 = lpv345.iter().filter_map(|p| p.sol_en).sum();
                let e2: f32 = lpv678.iter().filter_map(|p| p.sol_en).sum();
                let e1 = -e1 / n1;
                let e2 = -e2 / n2;
                n1w += if e1 > e2 { 1 } else { 0 };
                n2w += if e1 < e2 { 1 } else { 0 };
                n1e += e1;
                n2e += e2;
                println!("{s} - {f} t1:{e1:?} t2:{e2:?}");
            }
        }
    }
    println!("n1w:{n1w} n2w:{n2w} n1e:{n1e} n2e:{n2e}");
    Ok(())
}
