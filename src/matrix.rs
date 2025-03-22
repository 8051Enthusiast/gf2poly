use crate::Gf2Poly;

pub(crate) trait MatrixSubspace: core::ops::Mul<Self, Output = Self> + Sized {
    fn projection(matrix: Gf2Poly2x2Matrix) -> Self;
    fn quotient_matrix_projection(quotient: Gf2Poly) -> Self;
    fn identity() -> Self;
}

#[derive(Clone, Default)]
pub(crate) struct Gf2Poly2x2Matrix(pub Gf2Poly, pub Gf2Poly, pub Gf2Poly, pub Gf2Poly);

impl Gf2Poly2x2Matrix {
    pub fn multiply_simple(self, other: Gf2Poly2x2Matrix) -> Self {
        let a = &self.0 * &other.0 + &self.1 * &other.2;
        let b = &self.0 * &other.1 + &self.1 * &other.3;
        let c = &self.2 * &other.0 + &self.3 * &other.2;
        let d = &self.2 * &other.1 + &self.3 * &other.3;
        Self(a, b, c, d)
    }

    pub fn apply(&self, a: &Gf2Poly, b: &Gf2Poly) -> (Gf2Poly, Gf2Poly) {
        let x = &self.0 * a + &self.1 * b;
        let y = &self.2 * a + &self.3 * b;
        (x, y)
    }
}

impl core::ops::Mul<Gf2Poly2x2Matrix> for Gf2Poly2x2Matrix {
    type Output = Gf2Poly2x2Matrix;

    fn mul(self, other: Gf2Poly2x2Matrix) -> Self::Output {
        self.multiply_simple(other)
    }
}

impl core::fmt::Display for Gf2Poly2x2Matrix {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "[[{}, {}], [{}, {}]]", self.0, self.1, self.2, self.3)
    }
}

impl MatrixSubspace for Gf2Poly2x2Matrix {
    fn projection(matrix: Gf2Poly2x2Matrix) -> Self {
        matrix
    }

    fn quotient_matrix_projection(quotient: Gf2Poly) -> Self {
        Gf2Poly2x2Matrix(Gf2Poly::zero(), Gf2Poly::one(), Gf2Poly::one(), quotient)
    }

    fn identity() -> Self {
        Gf2Poly2x2Matrix(Gf2Poly::one(), Gf2Poly::zero(), Gf2Poly::zero(), Gf2Poly::one())
    }
}

/// This is so we can be generic in the hgcd routine for xgcd and regular gcd.
/// We don't want to calculate extra information for the gcd, so we just
/// throw it away by converting into this type.
pub(crate) struct TrivialSpace;

impl core::ops::Mul<Self> for TrivialSpace {
    type Output = Self;
    
    fn mul(self, _: Self) -> Self::Output {
        TrivialSpace
    }
}

impl MatrixSubspace for TrivialSpace {
    fn projection(_: Gf2Poly2x2Matrix) -> Self {
        Self
    }

    fn quotient_matrix_projection(_: Gf2Poly) -> Self {
        Self
    }

    fn identity() -> Self {
        Self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_multiply() {
        let a = Gf2Poly2x2Matrix(
            Gf2Poly::one(),
            Gf2Poly::x(),
            Gf2Poly::zero(),
            Gf2Poly::one(),
        );
        let b = Gf2Poly2x2Matrix(Gf2Poly::one(), Gf2Poly::x(), Gf2Poly::x(), Gf2Poly::x());
        assert_eq!((a.clone() * b.clone()).to_string(), "[[5, 6], [2, 2]]");
        assert_eq!((b * a).to_string(), "[[1, 0], [2, 6]]");
    }
}
