use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Solver {
    pub time_end: f32, // integration upper limit (starting from t=0)
    time_step: f32,    // time step Δt
    f_0: f32,          // initial value f(0)
    integral_f_0: f32, // initial value ∫f(0)dt
    mass: f32,         // mass
    spring: f32,       // spring constant
    amplitude: f32,    // amplitude of external force
    frequency: f32,    // frequency of external force
}

impl Solver {
    /*
     * We want:
     * - f(t)
     * - ∫f(t)dt for 0->t
     *
     * We don't know:
     * - exact function: f(t)
     * - exact function: ∫f(0)dt
     *
     * We know:
     * - derivetive of f(t): df(t)/dt
     * - initial value: f(0)
     * - initial value: ∫f(0)dt
     *
     * We have df(t)/f(t) and a way to approximate f(t+Δt) and ∫f(t+Δt)dt
     * - model: df(t)/dt = hoge_function(t, f(t), ∫f(t)dt)
     * - approx: f(t+Δt) = fuga_function(t, f(t), df(t)/dt)
     * - approx: ∫f(t+Δt)dt = piyo_function(t, ∫f(t)dt, f(t))
     *
     */
    pub fn new(
        time_end: f32,
        time_step: f32,
        f_0: f32,
        integral_f_0: f32,
        model_params: HashMap<String, f32>, // mass, spring, amplitude, frequency
    ) -> Solver {
        assert_eq!(
            model_params.len(),
            4,
            "model_params must be length 4: [mass, spring, amplitude, frequency]"
        );
        assert!(
            model_params.contains_key("mass"),
            "model_parameters must contain 'mass' key"
        );
        assert!(
            model_params.contains_key("spring"),
            "model_parameters must contain 'spring' key"
        );
        assert!(
            model_params.contains_key("amplitude"),
            "model_parameters must contain 'amplitude' key"
        );
        assert!(
            model_params.contains_key("frequency"),
            "model_parameters must contain 'frequency' key"
        );
        let mass = model_params["mass"];
        let spring = model_params["spring"];
        let amplitude = model_params["amplitude"];
        let frequency = model_params["frequency"];
        Solver {
            time_end,
            time_step,
            f_0,
            integral_f_0,
            mass,
            spring,
            amplitude,
            frequency,
        }
    }

    /* Analytical Solutions ******************************************************/
    pub fn function_t(&self, t: &f32) -> f32 {
        // Change here if analytical solution is known.
        // Analytical solution: f(t)
    }

    fn integral_function_t(&self, t: &f32) -> f32 {
        // Change here if analytical solution is known.
        // Analytical solution: ∫f(t)dt for 0 -> t
        // A = constant
        // omega^2 = spring_k / mass
        // sigma = constant
        // A cos(omega * t + sigma) + (ext_amp / omega^2) cos(ext_freq * t)
    }
    /* Analytical Solutions ********************************************************/

    /* Model ***********************************************************************/
    fn d_dt_function_t(
        &self,
        current_t: &f32,
        /*current_f_t: &f32,*/ current_integral: &f32,
    ) -> f32 {
        // Change this fuction to change the model.
        // df(t)/dt = -(spring_k/mass)*∫f(t)dt + ext_amp cos(ext_freq * t)
        -(current_integral * self.spring / self.mass)
            + self.amplitude * (self.frequency * current_t).cos()
    }
    /* Model ***********************************************************************/

    fn euler_next_step(
        &self,
        current_t: &f32,
        current_f_t: &f32,
        current_integral: &f32,
    ) -> (f32, f32) {
        // ∫f(t+Δt)dt = ∫f(t)dt + Δt * f(t)
        // f(t+Δt) = f(t) + Δt * df(t)/dt
        let next_integral = current_integral + self.time_step * current_f_t;
        let next_f_t =
            current_f_t + self.time_step * self.d_dt_function_t(current_t, current_integral);
        (next_f_t, next_integral)
    }

    fn euler_half_step(
        &self,
        current_t: &f32,
        current_f_t: &f32,
        current_integral: &f32,
    ) -> (f32, f32) {
        // ∫f(t+Δt)dt = ∫f(t)dt + (Δt/2) * f(t)
        // f(t+Δt) = f(t) + (Δt/2) * df(t)/dt
        let next_integral = current_integral + self.time_step * current_f_t / 2.0;
        let next_f_t =
            current_f_t + self.time_step * self.d_dt_function_t(current_t, current_integral) / 2.0;
        (next_f_t, next_integral)
    }

    // fn trapezoidal_next_step(
    //     &self, current_t: &f32, current_f_t: &f32, current_integral: &f32
    // ) -> f32 {
    //     // ∫f(t+Δt)dt = ∫f(t)dt + Δt/2 * [f(t) + f(t+Δt)]
    //     // f(t+Δt) = f(t) + Δt/2 * [df(t)/dt + df(t+Δt)/dt]
    //     // Using predictor-corrector approach
    // }

    fn rk4_next_step(
        &self,
        current_t: &f32,
        current_f_t: &f32,
        current_integral: &f32,
    ) -> (f32, f32) {
        // ∫f(t+Δt)dt = ∫f(t)dt + (Δt/6) * [k1_f(t) + 2*k2_f(t) + 2*k3_f(t) + k4_f(t)]
        // - k1_f(t) = func(t)
        // - k2_f(t) = func(t+Δt/2, ∫f(t)dt + k1_f(t) * Δt/2)
        // - k3_f(t) = func(t+Δt/2, ∫f(t)dt + k2_f(t) * Δt/2)
        // - k4_f(t) = func(t+Δt, ∫f(t)dt + k3_f(t) * Δt)
        // f(t+Δt) = f(t) + (Δt/6) * [k1_f'(t) + 2*k2_f'(t) + 2*k3_f'(t) + k4_f'(t)]
        // - k1_f'(t) = func(t, ∫f(t)dt)
        // - k2_f'(t) = func(t+Δt/2, ∫f(t)dt + k1_f'(t) * Δt/2)
        // - k3_f'(t) = func(t+Δt/2, ∫f(t)dt + k2_f'(t) * Δt/2)
        // - k4_f'(t) = func(t+Δt, ∫f(t)dt + k3_f'(t) * Δt)
        //
        let (euler_next_ft, _) = self.euler_next_step(current_t, current_f_t, current_integral);
        let (euler_half_ft, _) = self.euler_half_step(current_t, current_f_t, current_integral);
        let k_one_ft = current_t;
        let k_two_ft = euler_half_ft;
        let k_three_ft = euler_half_ft;
        let k_four_ft = euler_next_ft;
        let k_one_ddt_ft = self.d_dt_function_t(current_t, current_integral);
        let k_two_ddt_ft = self.d_dt_function_t(
            &(current_t + self.time_step / 2.0),
            &(current_integral + k_one_ddt_ft * self.time_step / 2.0),
        );
        let k_three_ddt_ft = self.d_dt_function_t(
            &(current_t + self.time_step / 2.0),
            &(current_integral + k_two_ddt_ft * self.time_step / 2.0),
        );
        let k_four_ddt_ft = self.d_dt_function_t(
            &(current_t + self.time_step),
            &(current_integral + k_three_ddt_ft * self.time_step),
        );

        let next_integral = current_integral
            + self.time_step * (k_one_ft + 2.0 * k_two_ft + 2.0 * k_three_ft + k_four_ft) / 6.0;
        let next_f_t = current_f_t
            + self.time_step
                * (k_one_ddt_ft + 2.0 * k_two_ddt_ft + 2.0 * k_three_ddt_ft + k_four_ddt_ft)
                / 6.0;
        (next_f_t, next_integral)
    }

    pub fn integrate_solution(&self, method: &str) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
        // Compute integral ∫f(t)dt from 0 to time_end using Euler method for ODE + trapezoidal rule for integration
        // method: &str
        // - "euler": euler method
        // // - "trap": trapezoid
        // - "rk4": Runge-Kutta 4
        // - "analytical": analytical
        let mut t: f32 = 0.0;
        let mut f_t: f32 = self.f_0;
        let mut integral: f32 = 0.0;

        let mut t_vec = vec![t];
        let mut f_t_vec = vec![f_t];
        let mut integral_vec = vec![integral];

        while t < self.time_end {
            // Compute f(t+Δt), ∫f(t+Δt)dt
            let f_t_next: f32;
            let integral_next: f32;
            if method == "euler" {
                (f_t_next, integral_next) = self.euler_next_step(&t, &f_t, &integral);
            // } else if method == "trap" {
            } else if method == "rk4" {
                (f_t_next, integral_next) = self.rk4_next_step(&t, &f_t, &integral);
            } else if method == "analytical" {
                f_t_next = self.function_t(&(t + self.time_step));
                integral_next = self.integral_function_t(&t);
            } else {
                panic!("Bad parameter; method: &str");
            }

            f_t = f_t_next;
            integral = integral_next;
            t += self.time_step;

            t_vec.push(t);
            f_t_vec.push(f_t);
            integral_vec.push(integral);
        }
        (t_vec, f_t_vec, integral_vec)
    }
}
