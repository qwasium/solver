use std::collections::HashMap;

use super::model::{Solver, SolverParam};

#[derive(Debug, Clone)]
pub struct Fall {
    pub condition: SolverParam,
    mass: f32,      // mass
    viscosity: f32, // viscosity coefficient
    gravity: f32,   // gravitational acceleration 9.8 m/s^2
}

impl Fall {
    pub fn new(
        time_end: f32,
        time_step: f32,
        init_acceleration: f32,
        init_velocity: f32,
        init_position: f32,
        model_params: HashMap<String, f32>, // {"mass": xx, "viscosity": yy}
    ) -> Self {
        assert_eq!(
            model_params.len(),
            2,
            "model_params must be length 2: [mass, viscosity]"
        );
        assert!(
            model_params.contains_key("mass"),
            "model_parameters must contain 'mass' key"
        );
        assert!(
            model_params.contains_key("viscosity"),
            "model_parameters must contain 'viscosity' key"
        );
        let mass = model_params["mass"];
        let viscosity = model_params["viscosity"];
        let condition = SolverParam {
            time_end,
            time_step,
            init_acceleration,
            init_velocity,
            init_position,
        };
        Self {
            condition,
            mass,
            viscosity,
            gravity: 9.8,
        }
    }
}

impl Solver for Fall {
    fn condition(&self) -> &SolverParam {
        &self.condition
    }
    // /* Analytical Solutions ******************************************************/
    // pub fn function_t(&self, t: &f32) -> f32 {
    //     // Change here if analytical solution is known.
    //     // Analytical solution: f(t)
    //     let mg_vis = self.mass * self.gravity / self.viscosity;
    //     mg_vis + (self.f_0 - mg_vis) * (-self.viscosity * t / self.mass).exp()
    // }
    //
    // fn integrate_function_t(&self, t: &f32) -> f32 {
    //     // Change here if analytical solution is known.
    //     // Analytical solution: ∫f(t)dt for 0 -> t
    //     let mg_vis = self.mass * self.gravity / self.viscosity;
    //     mg_vis * t
    //         - (self.f_0 - mg_vis)
    //             * (self.mass / self.viscosity)
    //             * ((-self.viscosity * t / self.mass).exp() - 1.0)
    // }
    // /* Analytical Solutions ********************************************************/
    /* Model ***********************************************************************/
    fn acceleration_t(&self, current_val: &f32, _current_position: &f32) -> f32 {
        // df(t)/dt = g - gamma*f(t)/mass
        self.gravity - self.viscosity * current_val / self.mass
    }
    /* Model ***********************************************************************/
}
