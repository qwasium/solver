use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Solver {
    pub time_end : f32, // integration upper limit (starting from t=0)
    time_step    : f32, // time step Δt
    f_0          : f32, // initial value f(0)
    mass         : f32, // mass
    viscosity    : f32, // viscosity coefficient
    gravity      : f32  // gravitational acceleration 9.8 m/s^2
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
    pub fn new(
        time_end: f32,
        time_step: f32,
        f_0: f32,
        model_params: HashMap<String, f32> // {"mass": xx, "viscosity": yy}
    ) -> Solver {
        assert_eq!(model_params.len(), 2, "model_params must be length 2: [mass, viscosity]");
        assert!(model_params.contains_key("mass"), "model_parameters must contain 'mass' key");
        assert!(model_params.contains_key("viscosity"), "model_parameters must contain 'viscosity' key");
        let mass = model_params["mass"];
        let viscosity = model_params["viscosity"];
        Solver {
            time_end, time_step, f_0,  mass, viscosity, gravity: 9.8
        }
    }

    /* Analytical Solutions ******************************************************/
    pub fn function_t(&self, t: &f32) -> f32 {
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
        let k_one = self.d_dt_function_t(current_val);
        let k_two = self.d_dt_function_t(&(current_val + self.time_step * k_one / 2.0));
        let k_three = self.d_dt_function_t(&(current_val + self.time_step * k_two / 2.0));
        let k_four = self.d_dt_function_t(&(current_val + self.time_step * k_three));
        current_val + (k_one + 2.0 * k_two + 2.0 * k_three + k_four) * self.time_step / 6.0
    }

    pub fn integrate_solution(&self, method: &str) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
        // Compute integral ∫f(t)dt from 0 to time_end using Euler method for ODE + trapezoidal rule for integration
        // method: &str
        // - "euler": euler method
        // - "trap": trapezoid
        // - "rk4": Runge-Kutta 4
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
