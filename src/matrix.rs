use crate::Gf2Poly;

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

    pub fn identity() -> Self {
        Gf2Poly2x2Matrix(
            Gf2Poly::one(),
            Gf2Poly::zero(),
            Gf2Poly::zero(),
            Gf2Poly::one(),
        )
    }

    pub fn apply(&self, a: &Gf2Poly, b: &Gf2Poly) -> (Gf2Poly, Gf2Poly) {
        let x = &self.0 * a + &self.1 * b;
        let y = &self.2 * a + &self.3 * b;
        (x, y)
    }
}

impl std::ops::Mul<Gf2Poly2x2Matrix> for Gf2Poly2x2Matrix {
    type Output = Gf2Poly2x2Matrix;

    fn mul(self, other: Gf2Poly2x2Matrix) -> Self::Output {
        self.multiply_simple(other)
    }
}

impl std::fmt::Display for Gf2Poly2x2Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[[{}, {}], [{}, {}]]", self.0, self.1, self.2, self.3)
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
