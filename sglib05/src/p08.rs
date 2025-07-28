use crate::p02::DrawLoadProf;
use crate::p02::LoadProf;
use crate::p02::SubLoadProf;
use std::collections::HashMap;
use std::error::Error;

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
