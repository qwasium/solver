
#[derive(Debug, Clone)]
struct Solver {
    time_end: f32,   // integration upper limit (starting from t=0)
    time_step: f32,
    f_0: f32,        // initial value f(0)
    mass: f32,       // mass
    viscosity: f32,  // viscosity coefficient
    gravity: f32     // gravitational acceleration 9.8 m/s^2
}

impl Solver {
    fn new(
        time_end: f32,
        time_step: f32,
        f_0: f32,
        mass: f32,
        viscosity: f32,
    ) -> Solver {
        Solver {
            time_end,
            time_step,
            f_0,
            mass,
            viscosity,
            gravity: 9.8
        }
    }
    
    fn function_t(&self, t: &f32) -> f32 {
        // Analytical solution: f(t)
        let mg_vis = self.mass * self.gravity / self.viscosity;
        mg_vis + (self.f_0 - mg_vis) * (-self.viscosity * t / self.mass).exp()
        // mg_vis * self.time_end - (self.f_0 - mg_vis) * mg_vis * ((-self.viscosity * self.time_end / self.mass).exp() - 1.0)
    }
    
    fn d_dt_function_t(&self, current_val: &f32) -> f32 {
        // df(t)/dt = g - gamma*f(t)/mass
        self.gravity - self.viscosity * current_val / self.mass
    }
    
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
    
    fn integrate_solution(&self, method: &str) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
        // Compute integral ∫f(t)dt from 0 to time_end using Euler method for ODE + trapezoidal rule for integration
        // method: &str
        // - "euler": euler method
        // - "trap": trapezoid
        let mut f_t: f32 = self.f_0;
        let mut area: f32 = 0.0;
        let mut t: f32 = 0.0;

        let mut t_vec = vec![t];
        let mut f_t_vec = vec![f_t];
        let mut area_vec = vec![area];
        
        while t < self.time_end {
            let f_t_next: f32;
            if method == "euler" {
                f_t_next = self.euler_next_step(&f_t);
            } else if method == "trap" {
                f_t_next = self.trapezoidal_next_step(&f_t);
            } else {
                panic!("Bad parameter; method: &str");
            }
            area += (f_t + f_t_next) * self.time_step / 2.0;

            f_t = f_t_next;
            t += self.time_step;

            t_vec.push(t);
            f_t_vec.push(f_t);
            area_vec.push(area);
        }
        (t_vec, f_t_vec, area_vec)
    }
    
    fn integrate_analytical(&self) -> f32 {
        // Compute integral of analytical solution for comparison
        let mut area = 0.0;
        let mut t = 0.0;

        while t < self.time_end {
            let f_t = self.function_t(&t);
            let f_t_next = self.function_t(&(t + self.time_step));
            area += (f_t + f_t_next) * self.time_step / 2.0;
            t += self.time_step;
        }
        area
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let conf = Solver::new(
        3.0,  // time_end
        0.05, // time_step
        2.0,  // f(0)
        1.0,  // mass
        1.0,  // viscosity
    );
    
    println!("{:?}", conf);
    println!();

    let (t_vec, euler_f_t, euler_area) = conf.integrate_solution("euler");
    let (_, trap_f_t, trap_area) = conf.integrate_solution("trap");
    
    // Compare final values f(time_end)
    println!("=== Final Values f({}) ===", conf.time_end);
    println!("Analytical solution: {:.6}", conf.function_t(&conf.time_end));
    println!("Euler method:        {:.6}", euler_f_t.last().unwrap());
    println!("Trapezoidal method:  {:.6}", trap_f_t.last().unwrap());
    println!();
    
    // Compare integrals ∫f(t)dt from 0 to time_end
    println!("=== Integrals ∫f(t)dt from 0 to {} ===", conf.time_end);
    println!("Analytical integral:           {:.6}", conf.integrate_analytical());
    println!("Euler ODE + Trap integration:  {:.6}", euler_area.last().unwrap());
    println!("Trap ODE + Trap integration:   {:.6}", trap_area.last().unwrap());
    println!();

    println!("t: {:?}", t_vec);
    Ok(())
}
