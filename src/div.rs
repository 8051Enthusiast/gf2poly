use crate::{BITS, Gf2Poly, Limb, limbs};
// the division algorithm here is described in
// http://people.csail.mit.edu/madhu/ST12/scribe/lect06.pdf

// naive algorithm for computing f^-1 mod x^8 for the lookup table
const fn const_byte_inverse(orig: u8) -> u8 {
    let mut product = orig;
    let mut result = 1;
    let mut i = 1;
    while i < 8 {
        if (product >> i) & 1 != 0 {
            product ^= orig << i;
            result |= 1 << i;
        }
        i += 1;
    }
    result
}

const BYTE_INVERSE: [Limb; 128] = {
    let mut res = [0; 128];
    let mut i = 0u8;
    while i < 128 {
        res[i as usize] = const_byte_inverse(2 * i + 1) as Limb;
        i += 1;
    }
    res
};

const BYTE_DEG: [u64; 128] = {
    let mut res = [0; 128];
    let mut i = 0;
    while i < 128 {
        res[i] = BITS - BYTE_INVERSE[i].leading_zeros() as u64 - 1;
        i += 1;
    }
    res
};

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
    let reduced = f.truncated(degree);
    let mut res = inverse_mod_power_impl(&reduced, degree);
    res.truncate_mut(degree);
    res
}

// it's assumed that `f.deg()` < `degree` and that f is invertible
fn inverse_mod_power_impl(f: &Gf2Poly, degree: u64) -> Gf2Poly {
    if degree <= 8 {
        let half = (f.limbs()[0] / 2) as usize;
        let inverse = BYTE_INVERSE[half];
        let deg = BYTE_DEG[half];
        return Gf2Poly {
            deg,
            limbs: limbs![inverse],
        };
    }

    let hdeg = degree.div_ceil(2);
    let lo = f.truncated(hdeg);
    let lo_inv = inverse_mod_power_impl(&lo, hdeg);
    // let half = x^hdeg, and let `lo` be f mod half, and lo_inv the inverse of lo mod half.
    // say hi_inv is such that inv = hi_inv * half + lo_inv is the inverse of f mod x^degree. then 
    // we equivalently have lo_inv = inv + hi_inv * half. let's calculate (modulo x^degree)
    // lo_inv * f = (inv + hi_inv * half) * f = 1 + half * hi_inv * f = 1 + half * hi_inv * lo (mod x^degree)
    // (note that we can turn the half * f into half * lo because the x^hdeg shifts the high bits of f out)
    // 
    // multiplying with lo_inv again we finally get:
    // lo_inv * (1 + half * hi_inv * lo) = lo_inv + half * lo_inv * lo * hi_inv = lo_inv + half * hi_inv = inv (mod x^degree)
    // 
    // also note that we get .square() basically for free in char 2.
    let mut inv = f * lo_inv.square();
    inv.truncate_mut(degree);
    inv
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
        (&self / rhs) * rhs + self
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
