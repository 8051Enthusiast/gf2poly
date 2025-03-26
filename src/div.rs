use crate::Gf2Poly;
// the division algorithm here is described in
// http://people.csail.mit.edu/madhu/ST12/scribe/lect06.pdf

// this calculates f^-1 mod x^degree, which is needed for the fast
// division algorithm.
fn inverse_mod_power(f: &Gf2Poly, degree: u64) -> Gf2Poly {
    // we do some preprocessing in this function so that we can
    // make assumptions about the input in inverse_mod_power_impl
    if degree == 0 {
        return Gf2Poly::zero();
    }
    if !f.eval(false) {
        panic!("f may not be divisible by x!");
    }
    let power_of_two_degree = degree.next_power_of_two();
    let reduced = f.truncated(power_of_two_degree);
    let mut res = inverse_mod_power_impl(&reduced, power_of_two_degree);
    res.truncate_mut(degree);
    res
}

// it's assumed that `f.deg()` < `degree` and degree is a power of two,
// also taht f is invertible
fn inverse_mod_power_impl(f: &Gf2Poly, degree: u64) -> Gf2Poly {
    if degree == 1 {
        return Gf2Poly::one();
    }

    // we split f into the halves lo and hi so that f = lo + hi * x^hdeg
    let hdeg = degree / 2;
    let lo = f.truncated(hdeg);
    let lo_inv = inverse_mod_power_impl(&lo, hdeg);
    // the product lo * lo_inv will be unit_hi * x^hdeg + 1 since the lower part has to
    // be lo * lo^-1 = 1 mod x^hdeg.
    let unit_hi = (&lo * &lo_inv) >> hdeg;
    let hi = f.clone() >> hdeg;
    let mut mul = &hi * &lo_inv;
    mul.truncate_mut(hdeg);
    mul += &unit_hi;
    let mut hi_inv = &lo_inv * &mul;
    hi_inv.truncate_mut(hdeg);
    (hi_inv << hdeg) + lo_inv
}

impl Gf2Poly {
    fn div_rev(&self, rhs: &Gf2Poly) -> Gf2Poly {
        // the trick here is that reversal of polynomial coefficients
        // is a multiplicative homomorphism (so rev(a) * rev(b) = rev(a * b))
        // the remainder `r` of a division f  = gq + r is at the lower part of the coefficients,
        // so if we reverse the polynomials we get
        // Rev(f) = Rev(q) * Rev(g) + x^(deg(f) - deg(r)) * Rev(r)
        // and if we calculate it modulo x^(deg(f) - deg(g)) we get to ignore the remainder
        // and can calculate the inverse modulo the power of x with the inverse_mod_power
        // function that uses hensel lifting above.
        let rev_lhs = self.reverse();
        let rev_rhs = rhs.reverse();
        let result_deg = self.deg() - rhs.deg();
        let rev_inv_rhs = inverse_mod_power(&rev_rhs, result_deg + 1);
        let mut rev_result = rev_lhs * rev_inv_rhs;
        rev_result.truncate_mut(result_deg + 1);
        let padding = result_deg - rev_result.deg();
        rev_result.reverse() << padding
    }

    fn divide(&self, rhs: &Gf2Poly) -> Gf2Poly {
        if rhs.is_zero() {
            panic!("Division by zero.");
        }
        if self.is_zero() || self.deg() < rhs.deg() {
            return Gf2Poly::zero();
        }
        if self.deg() == rhs.deg() {
            return Gf2Poly::one();
        }
        self.div_rev(rhs)
    }
}

impl core::ops::Div for &Gf2Poly {
    type Output = Gf2Poly;

    fn div(self, rhs: Self) -> Self::Output {
        self.divide(rhs)
    }
}

impl core::ops::Div for Gf2Poly {
    type Output = Gf2Poly;

    fn div(self, rhs: Self) -> Self::Output {
        &self / &rhs
    }
}

impl core::ops::Div<Gf2Poly> for &Gf2Poly {
    type Output = Gf2Poly;

    fn div(self, rhs: Gf2Poly) -> Self::Output {
        self / &rhs
    }
}

impl core::ops::Div<&Gf2Poly> for Gf2Poly {
    type Output = Gf2Poly;

    fn div(self, rhs: &Gf2Poly) -> Self::Output {
        &self / rhs
    }
}

impl core::ops::DivAssign<&Gf2Poly> for Gf2Poly {
    fn div_assign(&mut self, rhs: &Gf2Poly) {
        *self = &*self / rhs;
    }
}

impl core::ops::DivAssign<Gf2Poly> for Gf2Poly {
    fn div_assign(&mut self, rhs: Gf2Poly) {
        *self = &*self / &rhs;
    }
}

impl core::ops::Rem for &Gf2Poly {
    type Output = Gf2Poly;

    fn rem(self, rhs: Self) -> Self::Output {
        if rhs.is_zero() {
            return self.clone();
        }
        self + rhs * &(self / rhs)
    }
}

impl core::ops::Rem for Gf2Poly {
    type Output = Gf2Poly;

    fn rem(self, rhs: Self) -> Self::Output {
        if rhs.is_zero() {
            return self;
        }
        &self % &rhs
    }
}

impl core::ops::Rem<Gf2Poly> for &Gf2Poly {
    type Output = Gf2Poly;

    fn rem(self, rhs: Gf2Poly) -> Self::Output {
        self % &rhs
    }
}

impl core::ops::Rem<&Gf2Poly> for Gf2Poly {
    type Output = Gf2Poly;

    fn rem(self, rhs: &Gf2Poly) -> Self::Output {
        if rhs.is_zero() {
            return self;
        }
        (&self / rhs) * self + rhs
    }
}

impl core::ops::RemAssign<&Gf2Poly> for Gf2Poly {
    fn rem_assign(&mut self, rhs: &Gf2Poly) {
        if rhs.is_zero() {
            return;
        }
        *self = &*self % rhs;
    }
}

impl core::ops::RemAssign<Gf2Poly> for Gf2Poly {
    fn rem_assign(&mut self, rhs: Gf2Poly) {
        if rhs.is_zero() {
            return;
        }
        *self = &*self % &rhs;
    }
}

impl Gf2Poly {
    /// Calculates (quotient, remainder) of the division.
    pub fn divmod(&self, by: &Self) -> (Self, Self) {
        let quotient = self / by;
        let remainder = self + &quotient * by;
        (quotient, remainder)
    }

    /// Returns whether `self` divides `rhs`.
    pub fn divides(&self, rhs: &Self) -> bool {
        if rhs.is_zero() {
            return true;
        }
        (rhs % self).is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prop_assert_poly_eq;
    use proptest::prelude::*;

    #[test]
    fn zero_divisibility() {
        assert!(Gf2Poly::zero().divides(&Gf2Poly::zero()));
        assert!(Gf2Poly::one().divides(&Gf2Poly::zero()));
        assert!(!Gf2Poly::zero().divides(&Gf2Poly::one()));
    }

    proptest! {
        #[test]
        fn mod_inverse(mut a: Gf2Poly, degree in 1u64..128) {
            if !a.eval(false) {
                a += &Gf2Poly::one();
            }
            let a_inv = inverse_mod_power(&a, degree);
            prop_assert!(a_inv.deg() < degree);
            prop_assert!(a_inv.is_normalized());
            let one = (a * a_inv).truncated(degree);
            prop_assert_poly_eq!(one, Gf2Poly::one());
        }

        #[test]
        fn euclid_div(a: Gf2Poly, b: Gf2Poly) {
            prop_assume!(!b.is_zero());
            let q = &a / &b;
            let res = a + &q * &b;
            prop_assert!(q.is_normalized());
            prop_assert!(res.is_normalized());
            prop_assert!(res.is_zero() || res.deg() < b.deg());
        }
    }
}
