use crate::{Gf2Poly, matrix::Gf2Poly2x2Matrix};

fn hgcd(a0: &Gf2Poly, a1: &Gf2Poly) -> Gf2Poly2x2Matrix {
    eprintln!("{} {}", a0, a1);
    let halfdeg = a0.deg() / 2;
    if a1.deg() <= halfdeg {
        return Gf2Poly2x2Matrix::identity();
    }

    let a0_hi = a0.clone() >> halfdeg;
    let a1_hi = a1.clone() >> halfdeg;
    let first_matrix = hgcd(&a0_hi, &a1_hi);
    eprintln!("second half");
    drop((a0_hi, a1_hi));
    let (b0, b1) = first_matrix.apply(&a0, &a1);
    let (quotient, b2) = b0.divmod(&b1);
    drop(b0);
    let quotient_matrix =
        Gf2Poly2x2Matrix(Gf2Poly::zero(), Gf2Poly::one(), Gf2Poly::one(), quotient);
    let first_matrix_step = quotient_matrix * first_matrix;
    let quarterdeg = halfdeg / 2;
    let b1_hi = b1 >> quarterdeg;
    let b2_hi = b2 >> quarterdeg;
    let second_matrix = hgcd(&b1_hi, &b2_hi);
    drop((b1_hi, b2_hi));
    second_matrix * first_matrix_step
}

impl Gf2Poly {
    fn gcd_impl(&self, other: &Self) -> Self {
        if other.divides(&self) {
            return other.clone();
        }

        let (b0, b1) = hgcd(self, other).apply(self, other);
        let b2 = &b0 % &b1;
        if b2.is_zero() {
            return b1;
        }

        b1.gcd_impl(&b2)
    }

    pub fn gcd(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }

        let (a0, a1) = if other.deg() > self.deg() {
            (other, self)
        } else {
            (self, other)
        };

        a0.gcd_impl(a1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gcd1() {
        let a: Gf2Poly = "1caa4".parse().unwrap();
        let b: Gf2Poly = "18378".parse().unwrap();
        assert_eq!(a.gcd(&b).to_string(), "0");
    }
}
