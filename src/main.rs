
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
    
    fn solve_euler(&self) -> f32 {
        // Solve ODE to find f(time_end) using Euler's method
        let mut f_t = self.f_0;
        let mut t = 0.0;
        
        while t < self.time_end {
            f_t = self.euler_next_step(&f_t);
            t += self.time_step;
        }
        f_t
    }
    
    fn solve_trapezoidal(&self) -> f32 {
        // Solve ODE to find f(time_end) using trapezoidal method
        let mut f_t = self.f_0;
        let mut t = 0.0;
        
        while t < self.time_end {
            f_t = self.trapezoidal_next_step(&f_t);
            t += self.time_step;
        }
        f_t
    }
    
    fn integrate_euler_solution(&self) -> f32 {
        // Compute integral ∫f(t)dt from 0 to time_end using Euler method for ODE + trapezoidal rule for integration
        let mut f_t = self.f_0;
        let mut area = 0.0;
        let mut t = 0.0;
        
        while t < self.time_end {
            let f_t_next = self.euler_next_step(&f_t);
            area += (f_t + f_t_next) * self.time_step / 2.0;
            f_t = f_t_next;
            t += self.time_step;
        }
        area
    }
    
    fn integrate_trapezoidal_solution(&self) -> f32 {
        // Compute integral ∫f(t)dt from 0 to time_end using trapezoidal method for ODE + trapezoidal rule for integration
        let mut f_t = self.f_0;
        let mut area = 0.0;
        let mut t = 0.0;
        
        while t < self.time_end {
            let f_t_next = self.trapezoidal_next_step(&f_t);
            area += (f_t + f_t_next) * self.time_step / 2.0;
            f_t = f_t_next;
            t += self.time_step;
        }
        area
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
    
    // Compare final values f(time_end)
    println!("=== Final Values f({}) ===", conf.time_end);
    println!("Analytical solution: {:.6}", conf.function_t(&conf.time_end));
    println!("Euler method:        {:.6}", conf.solve_euler());
    println!("Trapezoidal method:  {:.6}", conf.solve_trapezoidal());
    println!();
    
    // Compare integrals ∫f(t)dt from 0 to time_end
    println!("=== Integrals ∫f(t)dt from 0 to {} ===", conf.time_end);
    println!("Analytical integral:           {:.6}", conf.integrate_analytical());
    println!("Euler ODE + Trap integration:  {:.6}", conf.integrate_euler_solution());
    println!("Trap ODE + Trap integration:   {:.6}", conf.integrate_trapezoidal_solution());
    
    Ok(())
}
