use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use csv::Writer;
use yaml_rust::YamlLoader;

use crate::model::Solver;
mod model;

fn read_yaml<P: AsRef<Path>>(fpath: P)
    -> Result<Vec<yaml_rust::Yaml>, Box<dyn Error>> {
    let f_str = fs::read_to_string(fpath)?;
    let docs = YamlLoader::load_from_str(&f_str)?;
    Ok(docs)
}

fn yaml_to_f32(yaml: &yaml_rust::Yaml) -> Option<f32> {
    // Convert yaml numeric to f32
    match yaml{
        yaml_rust::Yaml::Integer(i) => Some(*i as f32),
        yaml_rust::Yaml::Real(s) => s.parse().ok(),
        _ => None,
    }
}

fn yaml_to_f32_vec(yaml: &yaml_rust::Yaml) -> Option<Vec<f32>> {
    // Convert yaml array to Vec<f32>
    match yaml {
        yaml_rust::Yaml::Array(arr) => {
            let mut result = Vec::new();
            for item in arr {
                if let Some(val) = yaml_to_f32(item) {
                    result.push(val);
                } else {
                    return None; // Invalid item in array
                }
            }
            Some(result)
        }
        _ => None,
    }
}

fn write_csv<P: AsRef<Path>>(
    time_vec: &[f32], euler_vec: &[f32], trap_vec: &[f32], rk_vec: &[f32], anal_vec: &[f32],  filename: P
) -> Result<(), Box<dyn Error>> {
    assert_eq!(time_vec.len(), euler_vec.len(), "Input vector has unequal length.");
    assert_eq!(time_vec.len(), trap_vec.len(), "Input vector has unequal length.");
    assert_eq!(time_vec.len(), rk_vec.len(), "Input vector has unequal length.");
    assert_eq!(time_vec.len(), anal_vec.len(), "Input vector has unequal length.");

    if let Some(parent) = filename.as_ref().parent() {
        fs::create_dir_all(parent)?;
    }

    let mut wrt = Writer::from_path(filename)?;
    wrt.write_record(&["timestamp", "euler", "trapezoid", "rk4", "analytical"])?; // header
    for row in 0..time_vec.len() {
        wrt.write_record(
            &[
                time_vec[row].to_string(),  // col 0
                euler_vec[row].to_string(), // col 1
                trap_vec[row].to_string(),  // col 2
                rk_vec[row].to_string(),    // col 3
                anal_vec[row].to_string()   // col 4
            ]
        )?;
    }
    wrt.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // Read config/settings.yaml
    let conf_path = PathBuf::from("config").join("settings.yaml");
    let conf_docs = read_yaml(&conf_path)?;
    let doc = &conf_docs[0];
    // println!("{:?}", doc); // debug

    let conf = Solver::new(
        yaml_to_f32(&doc["time_end"]).expect("time_end must be numeric"),   // time_end
        yaml_to_f32(&doc["time_step"]).expect("time_step must be numeric"), // time_step
        yaml_to_f32(&doc["f_0"]).expect("f_0 must be numeric"),             // f(0)
        yaml_to_f32_vec(&doc["model_params"]).expect("model_params must be array of numerics")
    );

    println!("{:?}", conf);
    println!();

    let (t_vec, anal_f_t, anal_area) = conf.integrate_solution("analytical");
    let (_, euler_f_t, euler_area) = conf.integrate_solution("euler");
    let (_, trap_f_t, trap_area) = conf.integrate_solution("trap");
    let (_, rk_f_t, rk_area) = conf.integrate_solution("rk4");

    // Compare final values f(time_end)
    println!("=== Final Values f({}) ===", conf.time_end);
    println!("Analytical solution: {:.6}", conf.function_t(&conf.time_end));
    println!("Euler method:        {:.6}", euler_f_t.last().unwrap());
    println!("Trapezoidal method:  {:.6}", trap_f_t.last().unwrap());
    println!("Runge-Kutta method:  {:.6}", rk_f_t.last().unwrap());
    println!();

    // Compare integrals ∫f(t)dt from 0 to time_end
    println!("=== Integrals ∫f(t)dt from 0 to {} ===", conf.time_end);
    println!("Analytical integral:           {:.6}", anal_area.last().unwrap());
    println!("Euler ODE + Trap integration:  {:.6}", euler_area.last().unwrap());
    println!("Trap ODE + Trap integration:   {:.6}", trap_area.last().unwrap());
    println!("RK4 ODE + Trap integration:    {:.6}", rk_area.last().unwrap());
    println!();

    // Write result/*.csv
    let data_dir = PathBuf::from("result");
    let ft_path = data_dir.join(doc["ft_fname"].as_str().expect("ft_fname must be string"));
    let area_path = data_dir.join(doc["area_fname"].as_str().expect("area_fname must be string"));
    write_csv(&t_vec, &euler_f_t, &trap_f_t, &rk_f_t, &anal_f_t, &ft_path)?;
    write_csv(&t_vec, &euler_area, &trap_area, &rk_area, &anal_area, &area_path)?;
    Ok(())
}
