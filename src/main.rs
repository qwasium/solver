use csv::Writer;
use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use yaml_rust::YamlLoader;

mod model;
use model::fall::Fall;
use model::model::{Solver, SolverParam};
use model::vibration::Vibration;

#[derive(Debug)]
enum Model {
    Vibration(Vibration),
    Fall(Fall),
}

impl Model {
    fn integrate_solution(&self, method: &str) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
        match self {
            Model::Vibration(model) => model.integrate_solution(method),
            Model::Fall(model) => model.integrate_solution(method),
        }
    }
    fn condition(&self) -> &SolverParam {
        match self {
            Model::Vibration(model) => model.condition(),
            Model::Fall(model) => model.condition(),
        }
    }
}

fn read_yaml<P: AsRef<Path>>(fpath: P) -> Result<Vec<yaml_rust::Yaml>, Box<dyn Error>> {
    let f_str = fs::read_to_string(fpath)?;
    let docs = YamlLoader::load_from_str(&f_str)?;
    Ok(docs)
}

fn yaml_to_f32(yaml: &yaml_rust::Yaml) -> Option<f32> {
    // Convert yaml numeric to f32
    match yaml {
        yaml_rust::Yaml::Integer(i) => Some(*i as f32),
        yaml_rust::Yaml::Real(s) => s.parse().ok(),
        _ => None,
    }
}

fn yaml_to_f32_hashmap(yaml: &yaml_rust::Yaml) -> Option<HashMap<String, f32>> {
    // Convert yaml hash to HashMap<String, f32>
    match yaml {
        yaml_rust::Yaml::Hash(hash) => {
            let mut result = HashMap::new();
            for (key, value) in hash {
                if let Some(key_str) = key.as_str() {
                    if let Some(val) = yaml_to_f32(value) {
                        result.insert(key_str.to_string(), val);
                    } else {
                        return None; // Invalid value in hash
                    }
                } else {
                    return None; // Invalid key in hash
                }
            }
            Some(result)
        }
        _ => None,
    }
}

fn write_csv<P: AsRef<Path>>(
    time_vec: &[f32],
    euler_vec: &[f32],
    rk_vec: &[f32],
    // anal_vec: &[f32],
    filename: P,
) -> Result<(), Box<dyn Error>> {
    assert_eq!(
        time_vec.len(),
        euler_vec.len(),
        "Input vector has unequal length."
    );
    assert_eq!(
        time_vec.len(),
        rk_vec.len(),
        "Input vector has unequal length."
    );
    // assert_eq!(
    //     time_vec.len(),
    //     anal_vec.len(),
    //     "Input vector has unequal length."
    // );

    if let Some(parent) = filename.as_ref().parent() {
        fs::create_dir_all(parent)?;
    }

    let mut wrt = Writer::from_path(filename)?;
    wrt.write_record(&["timestamp", "euler", "rk4"])?; //, "analytical"])?; // header
    for row in 0..time_vec.len() {
        wrt.write_record(&[
            time_vec[row].to_string(),  // col 0
            euler_vec[row].to_string(), // col 1
            rk_vec[row].to_string(),    // col 3
                                        // anal_vec[row].to_string(),  // col 4
        ])?;
    }
    wrt.flush()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    // Read config
    let args: Vec<String> = env::args().collect();
    let conf_path = if args.len() > 1 {
        PathBuf::from(&args[1])
    } else {
        // default: config/settings.yaml
        PathBuf::from("config").join("settings.yaml")
    };
    let conf_docs = read_yaml(&conf_path)?;
    let doc = &conf_docs[0];
    // println!("{:?}", doc); // debug

    // simulation parameters
    let time_end = yaml_to_f32(&doc["time_end"]).expect("time_end must be numeric"); // time_end
    let time_step = yaml_to_f32(&doc["time_step"]).expect("time_step must be numeric"); // time_step
    let init_acc =
        yaml_to_f32(&doc["init_acceleration"]).expect("init_acceleration must be numeric"); // a(0)
    let init_vel = yaml_to_f32(&doc["init_velocity"]).expect("init_velocity must be numeric"); // v(0)
    let init_pos = yaml_to_f32(&doc["init_position"]).expect("init_position must be numeric"); // x(0)
    let params = yaml_to_f32_hashmap(&doc["model_params"])
        .expect("model_params must be hashmap of string and numerics");

    // load model
    let model = doc["model"].as_str().expect("model must be string");
    let conf = match model {
        "vibration" => Model::Vibration(Vibration::new(
            time_end, time_step, init_acc, init_vel, init_pos, params,
        )),
        "fall" => Model::Fall(Fall::new(
            time_end, time_step, init_acc, init_vel, init_pos, params,
        )),
        _ => {
            return Err(format!("Unknown model: {}", model).into());
        }
    };

    println!("{:?}", conf);
    println!();

    // let (t_vec, anal_acc, anal_vel, anal_pos) = conf.integrate_solution("analytical");
    let (t_vec, euler_acc, euler_vel, euler_pos) = conf.integrate_solution("euler");
    let (_, rk_acc, rk_vel, rk_pos) = conf.integrate_solution("rk4");

    // Compare final values f(time_end)
    println!("=== Final Values v({}) ===", conf.condition().time_end);
    // println!("Analytical solution: {:.6}", anal_vel.last().unwrap());
    println!("Euler method:        {:.6}", euler_vel.last().unwrap());
    println!("Runge-Kutta method:  {:.6}", rk_vel.last().unwrap());
    println!();

    // Compare integrals ∫f(t)dt from 0 to time_end
    println!(
        "=== Integrals x(t) = ∫v(t)dt from 0 to {} ===",
        conf.condition().time_end
    );
    // println!(
    //     "Analytical integral:           {:.6}",
    //     anal_pos.last().unwrap()
    // );
    println!(
        "Euler ODE + Trap integration:  {:.6}",
        euler_pos.last().unwrap()
    );
    println!(
        "RK4 ODE + Trap integration:    {:.6}",
        rk_pos.last().unwrap()
    );
    println!();

    // Write result/*.csv
    let data_dir = PathBuf::from("result");
    let acceleration_path =
        data_dir.join(doc["acc_fname"].as_str().expect("acc_fname must be string"));
    let velocity_path = data_dir.join(doc["vel_fname"].as_str().expect("vel_fname must be string"));
    let position_path = data_dir.join(doc["pos_fname"].as_str().expect("pos_fname must be string"));
    // write_csv(&t_vec, &euler_acc, &rk_acc, &anal_acc, &acceleration_path)?;
    // write_csv(&t_vec, &euler_acc, &rk_acc, &anal_acc, &velocity_path)?;
    // write_csv(&t_vec, &euler_pos, &rk_pos, &anal_pos, &position_path)?;
    write_csv(&t_vec, &euler_acc, &rk_acc, &acceleration_path)?;
    write_csv(&t_vec, &euler_acc, &rk_acc, &velocity_path)?;
    write_csv(&t_vec, &euler_pos, &rk_pos, &position_path)?;
    Ok(())
}
