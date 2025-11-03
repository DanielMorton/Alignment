use num_traits::Zero;

pub fn max<T: Zero + PartialOrd>(a: T, b: T) -> T {
    if a > b {
        a
    } else {
        b
    }
}

pub fn clamp_to_zero<T: Zero + PartialOrd>(score: T) -> T {
    max(score, T::zero())
}

pub trait Epsilon {
    fn epsilon() -> Self;
    fn fuzzy_equals(a: Self, b: Self) -> bool;
}

impl Epsilon for f64 {
    fn epsilon() -> Self {
        1e-6
    }

    fn fuzzy_equals(a: f64, b: f64) -> bool {
        (a - b).abs() < Self::epsilon()
    }
}

impl Epsilon for i32 {
    fn epsilon() -> Self {
        0
    }

    fn fuzzy_equals(a: i32, b: i32) -> bool {
        a == b
    }
}
