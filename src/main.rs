// use meval
// use saphyr::{Yaml, YamlLoader};
// use std::fs;

#[allow(dead_code)]
#[derive(Debug, Clone)]
struct Solver {
    time_start: f32,
    time_end: f32,
    time_step: f32,
    x_0: f32,
    // velocity_function: String,
}

impl Solver {

    fn new(time_start: f32, time_end: f32, time_step: f32, x_0: f32) -> Solver {
        Solver { time_start, time_end, time_step, x_0}
    }

    fn velocity_function (t: &f32) -> f32 {
        -t.exp()
    }

    fn trapezoidal_rule (&self) -> f32 {
        let mut area: f32 = 0.0;
        let mut t = self.time_start;
        while t < self.time_end {
            let v_t0 = Solver::velocity_function(&t);
            let v_t1 = Solver::velocity_function(&(t+self.time_step));  // TODO: cache this
            area += ((v_t0 + v_t1) * self.time_step) / 2.0;
            t += self.time_step;
        }
        area
    }

    // fn euler_method (&self) -> f32 {
    // }
    //
    // fn runge_kutta  (&self) -> f32 {
    // }
}

// fn print_type<T>(_: &T) { 
//     println!("{:?}", std::any::type_name::<T>());
// }



fn main() -> Result<(), Box<dyn std::error::Error>> {
    // let yaml_str = fs::read_to_string("config/settings.yaml")?;
    // println!("Raw YAML:\n{}", yaml_str);
    
    // let docs = Yaml::load_from_str(&yaml_str).unwrap();
    // let doc = &docs[0];

    // print_type(doc);
    
    // let conf = FunctionConf {
    //     t_start: doc["time_start"].as_f64().unwrap_or(0.0) as f32,
    //     t_end: doc["time_end"].as_f64().unwrap_or(0.0) as f32,
    //     t_step: doc["time_step"].as_f64().unwrap_or(0.1) as f32,
    //     x_0: doc["x_0"].as_f64().unwrap_or(0.0) as f32,
    //     v_func: doc["v_function"].as_str().unwrap_or("").to_string(),
    // };
    //
    // println!("Parsed config: {:?}", conf);
    // println!("Start time: {}", conf.t_start);
    // println!("End time: {}", conf.t_end);
    // println!("Time step: {}", conf.t_step);
    // println!("Initial x: {}", conf.x_0);
    // println!("Velocity function: {}", conf.v_func);

    let conf = Solver::new (0.0, 3.0, 0.1,0.0,
        // velocity_function: "exp(-t)".to_string(),
    );
    println!("{:?}", conf);
    // println!("{}", Solver::velocity_function(&conf.time_start));

    println!("Integration with trapezoidal rule: {}", conf.trapezoidal_rule());
    
    Ok(())
}
