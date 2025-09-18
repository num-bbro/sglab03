use std::env;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let now = std::time::SystemTime::now();
    let a1 = env::args().nth(1).unwrap_or("?".to_string());
    let a2 = env::args().nth(2).unwrap_or("?".to_string());
    let a3 = env::args().nth(3).unwrap_or("?".to_string());
    match a1.as_str() {
        "data_trans1" => sglib05::p08::data_trans1()?,
        "c01_chk_05" => sglib05::c04::c01_chk_05()?,
        "p08_calc_lp2" => sglib05::p08::p08_calc_lp2(a2.as_str(), a3.as_str())?,
        "p08_calc_lp1" => sglib05::p08::p08_calc_lp1(a2.as_str())?,
        "p08_draw_01" => sglib05::p08::p08_draw_01(a2.as_str(), a3.as_str())?,
        "p08_draw_slp" => sglib05::p08::p08_draw_slp(a2.as_str(), a3.as_str())?,
        "p07_lp_pro" => sglib05::p07::p07_lp_pro(a2.as_str())?,
        "p03_calc_lp3" => sglib05::p03::p03_calc_lp3(a2.as_str())?,
        "p03_calc_lp2" => sglib05::p03::p03_calc_lp2(a2.as_str(), a3.as_str())?,
        "c04_chk_lp_01" => sglib05::c04::c04_chk_lp_01()?,
        "c01_chk_04" => sglib05::c04::c01_chk_04()?,
        "c01_chk_03" => sglib05::c04::c01_chk_03()?,
        "c01_chk_02" => sglib05::c04::c01_chk_02()?,
        "c01_chk_01" => sglib05::c04::c01_chk_01()?,
        "p05_ana_1" => sglib05::p05::p05_ana_1()?,
        "p06_amp_chk2" => sglib05::p06::p06_amp_chk2()?,
        "p01_ana_test3" => sglib05::p01::p01_ana_test3()?,
        "p01_ana_test2" => sglib05::p01::p01_ana_test2()?,
        "p02_read_lp" => sglib05::p02::p02_read_lp(a2.as_str())?,
        "p02_lp_pro" => sglib05::p02::p02_lp_pro(a2.as_str())?,
        "p02_draw_lp" => sglib05::p02::p02_draw_lp(a2.as_str())?,
        "p03_calc_lp" => sglib05::p03::p03_calc_lp(a2.as_str())?,
        "p03_draw_slp" => sglib05::p03::p03_draw_slp(a2.as_str(), a3.as_str())?,
        "p03_draw_all" => sglib05::p03::p03_draw_all()?,
        "p04_form_sub" => sglib05::p04::p04_form_sub()?,
        n => {
            println!("'{}' NG command", n);
        }
    }
    let se = now.elapsed().unwrap().as_secs();
    let mi = se / 60;
    println!("time {se} sec = {mi} min");
    Ok(())
}
