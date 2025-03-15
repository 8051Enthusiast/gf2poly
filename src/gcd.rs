use crate::{Gf2Poly, matrix::Gf2Poly2x2Matrix};

fn iter_hgcd(a0: &mut Gf2Poly, a1: &mut Gf2Poly) -> Gf2Poly2x2Matrix {
    let halfdeg = a0.deg() / 2;
    debug_assert!(a0.deg() >= a1.deg());

    let mut a00 = Gf2Poly::one();
    let mut a01 = Gf2Poly::zero();
    let mut a10 = Gf2Poly::zero();
    let mut a11 = Gf2Poly::one();

    loop {
        if a1.deg() <= halfdeg {
            break Gf2Poly2x2Matrix(a00, a01, a10, a11);
        }

        let (q, r) = a0.divmod(&a1);
        std::mem::swap(a0, a1);
        *a1 = r;

        std::mem::swap(&mut a00, &mut a10);
        std::mem::swap(&mut a01, &mut a11);

        a00 += &(&a10 * &q);
        a01 += &(&a11 * &q);
    }
}

fn hgcd(a0: &mut Gf2Poly, a1: &mut Gf2Poly) -> Gf2Poly2x2Matrix {
    let halfdeg = a0.deg() / 2;
    debug_assert!(a0.deg() >= a1.deg());
    if a1.deg() <= halfdeg {
        return Gf2Poly2x2Matrix::identity();
    }

    if a0.deg() < 4 {
        return iter_hgcd(a0, a1);
    }

    let mut a0_hi = a0.clone() >> halfdeg;
    let mut a1_hi = a1.clone() >> halfdeg;
    let first_matrix = hgcd(&mut a0_hi, &mut a1_hi);
    drop((a0_hi, a1_hi));
    let (b0, b1) = first_matrix.apply(&a0, &a1);
    if b1.deg() <= halfdeg {
        *a0 = b0;
        *a1 = b1;
        return first_matrix;
    }
    let (quotient, b2) = b0.divmod(&b1);
    drop(b0);
    let quotient_matrix =
        Gf2Poly2x2Matrix(Gf2Poly::zero(), Gf2Poly::one(), Gf2Poly::one(), quotient);
    let first_matrix_step = quotient_matrix * first_matrix;
    let quarterdeg = halfdeg / 2;
    let mut b1_hi = b1.clone() >> quarterdeg;
    let mut b2_hi = b2.clone() >> quarterdeg;
    let second_matrix = hgcd(&mut b1_hi, &mut b2_hi);
    drop((b1_hi, b2_hi));
    let (res_a0, res_a1) = second_matrix.apply(&b1, &b2);
    *a0 = res_a0;
    *a1 = res_a1;
    second_matrix * first_matrix_step
}

impl Gf2Poly {
    fn gcd_impl(&mut self, other: &mut Self) -> Self {
        let a = self;
        let b = other;

        while !b.is_zero() {
            debug_assert!(a.deg() >= b.deg());
            if b.divides(a) {
                return std::mem::take(b);
            }

            hgcd(a, b);
            if b.is_zero() {
                return std::mem::take(a);
            }
            let c = &*a % &*b;

            *a = std::mem::take(b);
            *b = c;
        }

        std::mem::take(a)
    }

    pub fn gcd(self, other: Self) -> Self {
        if self.is_zero() {
            return other;
        }
        if other.is_zero() {
            return self;
        }

        let (mut a0, mut a1) = if other.deg() > self.deg() {
            (other, self)
        } else {
            (self, other)
        };

        a0.gcd_impl(&mut a1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn gcd1() {
        let a: Gf2Poly = "1caa4".parse().unwrap();
        let b: Gf2Poly = "18378".parse().unwrap();
        assert_eq!(a.gcd(b).to_string(), "39c");
    }

    #[test]
    fn gcd2() {
        let a: Gf2Poly = "3eb64ce7d0c7fdce57504c4bf289023c4f7d42408b0d743ea5c8afb6b745f40c".parse().unwrap();
        let b: Gf2Poly = "4a3f022600075687ff0a08d032b36f22a7dff447673d5b6d4c9409be3ecd31d4".parse().unwrap();
        assert_eq!(a.gcd(b).to_string(), "14be7a0d84bbf65776cec6e7ebbb4c50c");
    }

    proptest! {
        #[test]
        fn gcd(a: Gf2Poly, b: Gf2Poly, c: Gf2Poly) {
            let amul = &a * &c;
            let bmul = &b * &c;
            let gcd = amul.clone().gcd(bmul.clone());
            prop_assert!(gcd.divides(&amul));
            prop_assert!(gcd.divides(&bmul));
            prop_assert!(c.divides(&gcd));
        }
    }
}
