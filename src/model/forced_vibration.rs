// use num::pow;
use std::collections::HashMap;

use super::model::{Solver, SolverParam};

#[derive(Debug, Clone)]
pub struct ForcedVibration {
    pub condition: SolverParam,
    mass: f32,      // mass
    spring: f32,    // spring constant
    amplitude: f32, // amplitude of external force
    frequency: f32, // frequency of external force
}

impl ForcedVibration {
    pub fn new(
        time_end: f32,
        time_step: f32,
        init_acceleration: f32,
        init_velocity: f32,
        init_position: f32,
        model_params: HashMap<String, f32>, // mass, spring, amplitude, frequency
    ) -> Self {
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
            spring,
            amplitude,
            frequency,
        }
    }
}

impl Solver for ForcedVibration {
    fn condition(&self) -> &SolverParam {
        &self.condition
    }

    // /* Analytical Solutions ******************************************************/
    // fn analytical_t(&self, t: &f32) -> (f32, f32, f32) {
    //     // A = ext_amp / (spring_k - mass * ext_freq^2)
    //     // delta = 0
    //     // x(t) = A cos(sqrt(spring_k / mass) * t + delta) + (ext_amp / (spring_k - mass * ext_freq^2)) * cos(ext_freq * t)
    //     let const_a = self.amplitude / (self.spring - self.mass * pow(self.frequency, 2));
    //     let delta: f32 = 0.0;
    //     let pos = const_a * ((self.spring / self.mass).sqrt() * t + delta).cos()
    //         + const_a * (self.frequency * t).cos();
    //     let vel = -const_a
    //         * (self.spring / self.mass).sqrt()
    //         * ((self.spring / self.mass).sqrt() * t + delta).sin()
    //         - const_a * self.frequency * (self.frequency * t).sin();
    //     let acc = -const_a
    //         * (self.spring / self.mass)
    //         * ((self.spring / self.mass).sqrt() * t + delta).cos()
    //         - const_a * pow(self.frequency, 2) * (self.frequency * t).sin();
    //     (acc, vel, pos)
    // }
    // /* Analytical Solutions ********************************************************/
    /* Model ***********************************************************************/
    fn acceleration_t(&self, current_t: &f32, current_position: &f32) -> f32 {
        // a(t)
        //   = -(spring_k / mass) * x(t)dt + (ext_amp / mass) * cos(ext_freq * t)
        //   = -(spring_k / mass) * ∫v(t)dt + (ext_amp / mass) * cos(ext_freq * t)
        -(current_position * self.spring / self.mass)
            + (self.amplitude / self.mass) * (self.frequency * current_t).cos()
    }
    /* Model ***********************************************************************/
}
