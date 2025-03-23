use crate::{
    Gf2Poly,
    matrix::{Gf2Poly2x2Matrix, MatrixSubspace, TrivialSpace},
};

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
        core::mem::swap(a0, a1);
        *a1 = r;

        core::mem::swap(&mut a00, &mut a10);
        core::mem::swap(&mut a01, &mut a11);

        a10 += &(&a00 * &q);
        a11 += &(&a01 * &q);
    }
}

fn hgcd<M: MatrixSubspace>(a0: &mut Gf2Poly, a1: &mut Gf2Poly) -> M {
    let halfdeg = a0.deg() / 2;
    debug_assert!(a0.deg() >= a1.deg());
    if a1.deg() <= halfdeg {
        return M::identity();
    }

    if a0.deg() < 4 {
        return M::projection(iter_hgcd(a0, a1));
    }

    let mut a0_hi = a0.clone() >> halfdeg;
    let mut a1_hi = a1.clone() >> halfdeg;
    let first_matrix = hgcd::<Gf2Poly2x2Matrix>(&mut a0_hi, &mut a1_hi);
    drop((a0_hi, a1_hi));
    let (b0, b1) = first_matrix.apply(&a0, &a1);
    if b1.deg() <= halfdeg {
        *a0 = b0;
        *a1 = b1;
        return M::projection(first_matrix);
    }
    let (quotient, b2) = b0.divmod(&b1);
    drop(b0);
    let quotient_matrix = M::quotient_matrix_projection(quotient);
    let first_matrix_step = quotient_matrix * M::projection(first_matrix);
    let quarterdeg = halfdeg / 2;
    let mut b1_hi = b1.clone() >> quarterdeg;
    let mut b2_hi = b2.clone() >> quarterdeg;
    let second_matrix = hgcd::<Gf2Poly2x2Matrix>(&mut b1_hi, &mut b2_hi);
    drop((b1_hi, b2_hi));
    let (res_a0, res_a1) = second_matrix.apply(&b1, &b2);
    let res = M::projection(second_matrix) * first_matrix_step;
    *a0 = res_a0;
    *a1 = res_a1;
    res
}

impl Gf2Poly {
    fn gcd_impl<M: MatrixSubspace>(&mut self, other: &mut Self) -> (Self, M) {
        let a = self;
        let b = other;
        let mut acc = M::identity();

        while !b.is_zero() {
            debug_assert!(a.deg() >= b.deg());
            let (quotient, remainder) = a.divmod(b);
            if remainder.is_zero() {
                acc = M::quotient_matrix_projection(quotient) * acc;
                return (core::mem::take(b), acc);
            }

            acc = hgcd::<M>(a, b) * acc;
            if b.is_zero() {
                return (core::mem::take(a), acc);
            }
            let (quotient, remainder) = a.divmod(b);
            acc = M::quotient_matrix_projection(quotient) * acc;

            *a = core::mem::take(b);
            *b = remainder;
        }

        (core::mem::take(a), acc)
    }

    /// Calculates the greatest common divisor of `self` and `other`.
    /// The result has the property that it divides `self` and `other` and
    /// it divides no other element that has that property.
    ///
    /// # Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a: Gf2Poly = "1caa4".parse().unwrap();
    /// let b: Gf2Poly = "18378".parse().unwrap();
    /// assert_eq!(a.gcd(b).to_string(), "39c");
    /// ```
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

        a0.gcd_impl::<TrivialSpace>(&mut a1).0
    }

    /// Calculates the extended greatest common divisor of `self` and `other`.
    /// Besides returning what `gcd` would return, it also returns `[x, y]`
    /// such that `self * x + other * y == gcd`.
    ///
    /// # Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let a: Gf2Poly = "1caa4".parse().unwrap();
    /// let b: Gf2Poly = "18378".parse().unwrap();
    /// let (gcd, [x, y]) = a.clone().xgcd(b.clone());
    /// assert_eq!(gcd.to_string(), "39c");
    /// assert_eq!((a * x + b * y).to_string(), "39c");
    /// ```
    pub fn xgcd(self, other: Self) -> (Gf2Poly, [Gf2Poly; 2]) {
        match (self.is_zero(), other.is_zero()) {
            (true, true) => return (Gf2Poly::zero(), [Gf2Poly::zero(), Gf2Poly::zero()]),
            (true, false) => return (other.clone(), [Gf2Poly::zero(), Gf2Poly::one()]),
            (false, true) => return (self.clone(), [Gf2Poly::one(), Gf2Poly::zero()]),
            _ => {}
        }

        let swap = other.deg() > self.deg();
        let (mut a0, mut a1) = if swap { (other, self) } else { (self, other) };

        let (gcd, matrix) = a0.gcd_impl::<Gf2Poly2x2Matrix>(&mut a1);
        let mut x = matrix.0;
        let mut y = matrix.1;
        if swap {
            core::mem::swap(&mut x, &mut y);
        }
        (gcd, [x, y])
    }
}

#[cfg(test)]
mod tests {
    use crate::prop_assert_poly_eq;

    use super::*;
    use proptest::prelude::*;

    #[test]
    fn gcd_big() {
        let a: Gf2Poly = "3eb64ce7d0c7fdce57504c4bf289023c4f7d42408b0d743ea5c8afb6b745f40c"
            .parse()
            .unwrap();
        let b: Gf2Poly = "4a3f022600075687ff0a08d032b36f22a7dff447673d5b6d4c9409be3ecd31d4"
            .parse()
            .unwrap();
        assert_eq!(a.gcd(b).to_string(), "14be7a0d84bbf65776cec6e7ebbb4c50c");
    }

    #[test]
    fn xgcd_big() {
        let a: Gf2Poly = "f2341b2123ad24c3b2e829e".parse().unwrap();
        let b: Gf2Poly = "1052649a".parse().unwrap();
        let (gcd, [x, y]) = a.clone().xgcd(b.clone());
        assert_eq!(a * x + b * y, gcd);
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

        #[test]
        fn xgcd(a: Gf2Poly, b: Gf2Poly, c: Gf2Poly) {
            let amul = &a * &c;
            let bmul = &b * &c;
            let (gcd, [x, y]) = amul.clone().xgcd(bmul.clone());
            prop_assert_poly_eq!(amul * x + bmul * y, gcd);
        }
    }
}
