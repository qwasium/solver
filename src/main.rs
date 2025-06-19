use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use csv::Writer;
use yaml_rust::YamlLoader;

#[derive(Debug, Clone)]
struct Solver {
    time_end : f32, // integration upper limit (starting from t=0)
    time_step: f32, // time step Δt
    f_0      : f32, // initial value f(0)
    mass     : f32, // mass
    viscosity: f32, // viscosity coefficient
    gravity  : f32  // gravitational acceleration 9.8 m/s^2
}

impl Solver {
    /*
     * We want:
     * - f(t)
     * - ∫f(t)dt for 0->t
     *
     * We don't know:
     * - exact function: f(t)
     *
     * We know:
     * - derivetive of f(t): df(t)/dt
     * - initial value: f(0)
     *
     * We have a model of f(current) & f'(current) -> f(next)
     * - model: f(t+Δt) = some_function(f(t), f'(t))
     *
     */
    fn new(
        time_end: f32,
        time_step: f32,
        f_0: f32,
        mass: f32,
        viscosity: f32,
    ) -> Solver {
        Solver {time_end, time_step, f_0, mass, viscosity, gravity: 9.8}
    }

    /* Analytical Solutions ******************************************************/
    fn function_t(&self, t: &f32) -> f32 {
        // Change here if analytical solution is known.
        // Analytical solution: f(t)
        let mg_vis = self.mass * self.gravity / self.viscosity;
        mg_vis + (self.f_0 - mg_vis) * (-self.viscosity * t / self.mass).exp()
    }

    fn integrate_function_t(&self, t: &f32) -> f32 {
        // Change here if analytical solution is known.
        // Analytical solution: ∫f(t)dt for 0 -> t
        let mg_vis = self.mass * self.gravity / self.viscosity;
        mg_vis * t - (self.f_0 - mg_vis) * (self.mass / self.viscosity) * ((-self.viscosity * t / self.mass).exp() - 1.0)
    }
    /* Analytical Solutions ********************************************************/

    /* Model ***********************************************************************/
    fn d_dt_function_t(&self, current_val: &f32) -> f32 {
        // Change this fuction to change the model.
        // df(t)/dt = g - gamma*f(t)/mass
        self.gravity - self.viscosity * current_val / self.mass
    }
    /* Model ***********************************************************************/

    fn euler_next_step(&self, current_val: &f32) -> f32 {
        // f(t+Δt) = f(t) + Δt * df(t)/dt
        current_val + self.time_step * self.d_dt_function_t(current_val)
    }

    fn trapezoidal_next_step(&self, current_val: &f32) -> f32 {
        // f(t+Δt) = f(t) + Δt/2 * [df(t)/dt + df(t+Δt)/dt]
        // Using predictor-corrector approach
        let predictor = current_val + self.time_step * self.d_dt_function_t(current_val);
        current_val + self.time_step * (self.d_dt_function_t(current_val) + self.d_dt_function_t(&predictor)) / 2.0
    }

    fn rk4_next_step(&self, current_val: &f32) -> f32 {
        // f(t+Δt) = f(t) + (Δt/6) * (k_1 + 2*k_2 + 2*k_3 + k_4)
        // k_1: df(t)/dt
        // k_2: 
        let k_one = self.d_dt_function_t(current_val);
        let k_two = self.d_dt_function_t(&(current_val + self.time_step * k_one / 2.0));
        let k_three = self.d_dt_function_t(&(current_val + self.time_step * k_two / 2.0));
        let k_four = self.d_dt_function_t(&(current_val + self.time_step * k_three));
        current_val + (k_one + 2.0 * k_two + 2.0 * k_three + k_four) * self.time_step / 6.0
    }

    fn integrate_solution(&self, method: &str) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
        // Compute integral ∫f(t)dt from 0 to time_end using Euler method for ODE + trapezoidal rule for integration
        // method: &str
        // - "euler": euler method
        // - "trap": trapezoid
        // - "analytical": analytical 
        let mut t: f32 = 0.0;
        let mut f_t: f32 = self.f_0;
        let mut area: f32 = 0.0;

        let mut t_vec = vec![t];
        let mut f_t_vec = vec![f_t];
        let mut area_vec = vec![area];

        while t < self.time_end {
            // Compute f(t+Δt)
            let f_t_next: f32;
            if method == "euler" {
                f_t_next = self.euler_next_step(&f_t);
                area += (f_t + f_t_next) * self.time_step / 2.0;
            } else if method == "trap" {
                f_t_next = self.trapezoidal_next_step(&f_t);
                area += (f_t + f_t_next) * self.time_step / 2.0;
            } else if method == "rk4" {
                f_t_next = self.rk4_next_step(&f_t);
                area += (f_t + f_t_next) * self.time_step / 2.0;
            } else if method == "analytical" {
                f_t_next = self.function_t(&(t + self.time_step));
                area = self.integrate_function_t(&t);
            } else {
                panic!("Bad parameter; method: &str");
            }

            f_t = f_t_next;
            t += self.time_step;

            t_vec.push(t);
            f_t_vec.push(f_t);
            area_vec.push(area);
        }
        (t_vec, f_t_vec, area_vec)
    }
}

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
        yaml_to_f32(&doc["mass"]).expect("mass must be numeric"),           // mass
        yaml_to_f32(&doc["viscosity"]).expect("viscosity must be numeric"), // viscosity
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
