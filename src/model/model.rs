#[derive(Debug, Clone)]
pub struct SolverParam {
    pub time_end: f32,          // integration upper limit (starting from t=0)
    pub time_step: f32,         // time step Δt
    pub init_acceleration: f32, // dv(0)/dt = a(0), t=0
    pub init_velocity: f32,     // v(0), t=0
    pub init_position: f32,     // ∫v(0)dt = x(0), t=0
}

pub trait Solver {
    /*
     * - acceleration: a(t) = dv(t)/dt
     * - velocity    : v(t)
     * - postion     : x(t) = ∫v(t)dt
     *
     * We want to compute:
     * - velocity at given t = v(t)
     * - postion  at given t = x(t) = ∫v(t)dt for 0->t
     * Where exact function of v(t) and x(t) are unknown.
     *
     * We know:
     * - initial acceleration : a(0)
     * - initial velocity     : v(0)
     * - initial postion      : x(0) = ∫v(0)dt
     * - model of acceleration: a(t) = some_function(t, v(t))
     *
     * We have df(t)/f(t) and a way to approximate v(t+Δt) and x(t+Δt)
     * - approximation: a(t+Δt) = some_function(t, v(t))
     *   - euler method
     *   - runge kutta method
     * - v(t+Δt) = v(t) + Δt * a(t)
     * - x(t+Δt) = x(t) + Δt * v(t)
     *
     */

    fn condition(&self) -> &SolverParam;

    // fn analytical_t(&self, t: &f32) -> (f32, f32, f32); // -> (acc, vel, pos)
    // Change here if analytical solution is known.
    // Analytical solution:
    // - x(t) = ∫v(t)dt for 0 -> t
    // - v(t) = dx(t)/dt
    // - a(t) = dv(t)/dt

    fn acceleration_t(&self, current_t: &f32, current_position: &f32) -> f32;
    // Change this fuction to change the model.
    // a(t) = dv(t)/dt = f(t, v(t))

    fn euler_next_step(
        &self,
        current_t: &f32,
        current_acceleration: &f32,
        current_velocity: &f32,
        current_position: &f32,
    ) -> (f32, f32, f32) {
        // v(t+Δt) = v(t) + Δt * f(t, a(t))
        // x(t+Δt) = x(t) + Δt * F(t, v(t))
        let cond = self.condition();
        let next_acc = current_acceleration
            + cond.time_step * self.acceleration_t(current_t, current_position);
        let next_vel = current_velocity + cond.time_step * (*current_acceleration + next_acc) / 2.0;
        let next_pos = current_position + cond.time_step * (current_velocity + next_vel) / 2.0;
        (next_acc, next_vel, next_pos)
    }

    fn rk4_next_step(
        &self,
        current_t: &f32,
        current_acceleration: &f32,
        current_velocity: &f32,
        current_position: &f32,
    ) -> (f32, f32, f32) {
        // runge-kutta and trapezoidal
        // a(t+Δt)dt = a(t) + (Δt/6) * [k1 + 2 * k2 + 2 * k3 + k4]
        // - k1 = f(t, v(t))
        // - k2 = f(t+Δt/2, v(t) + k1 * Δt/2)
        // - k3 = f(t+Δt/2, v(t) + k2 * Δt/2)
        // - k4 = f(t+Δt, v(t) + k3 * Δt)
        //
        let cond = self.condition();
        let half_timestep = cond.time_step / 2.0;
        let k_one = self.acceleration_t(current_t, current_position);
        let k_two = self.acceleration_t(
            &(current_t + half_timestep),
            &(current_position + &half_timestep * &k_one * (current_t + half_timestep)),
        );
        let k_three = self.acceleration_t(
            &(current_t + half_timestep),
            &(current_position + &half_timestep * &k_two * (current_t + half_timestep)),
        );
        let k_four = self.acceleration_t(
            &(current_t + cond.time_step),
            &(current_position + cond.time_step * &k_three * (current_t + cond.time_step)),
        );
        let next_acc = current_acceleration
            + (cond.time_step / 6.0) * (k_one + 2.0 * k_two + 2.0 * k_three + k_four);

        let next_vel = current_velocity + cond.time_step * (current_acceleration + next_acc) / 2.0;
        let next_pos = current_position + cond.time_step * (current_velocity + next_vel) / 2.0;
        (next_acc, next_vel, next_pos)
    }

    fn integrate_solution(&self, method: &str) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
        // Compute integral ∫f(t)dt from 0 to time_end using Euler method for ODE + trapezoidal rule for integration
        // method: &str
        // - "euler": euler method
        // // - "trap": trapezoid
        // - "rk4": Runge-Kutta 4
        // // - "analytical": analytical
        let cond = self.condition();
        let mut t: f32 = 0.0;
        let mut acc: f32 = cond.init_acceleration;
        let mut vel: f32 = cond.init_velocity;
        let mut pos: f32 = cond.init_position;

        let mut t_vec = vec![t];
        let mut acc_vec = vec![acc];
        let mut vel_vec = vec![vel];
        let mut pos_vec = vec![pos];

        while t < cond.time_end {
            // Compute f(t+Δt), ∫f(t+Δt)dt
            let acc_next: f32;
            let vel_next: f32;
            let pos_next: f32;
            if method == "euler" {
                (acc_next, vel_next, pos_next) = self.euler_next_step(&t, &acc, &vel, &pos);
            } else if method == "rk4" {
                (acc_next, vel_next, pos_next) = self.rk4_next_step(&t, &acc, &vel, &pos);
            // } else if method == "analytical" {
            //     (acc_next, vel_next, pos_next) = self.analytical_t(&t);
            } else {
                panic!("Bad parameter; method: &str");
            }

            acc = acc_next;
            vel = vel_next;
            pos = pos_next;
            t += cond.time_step;

            t_vec.push(t);
            acc_vec.push(acc);
            vel_vec.push(vel);
            pos_vec.push(pos);
        }
        (t_vec, acc_vec, vel_vec, pos_vec)
    }
}
