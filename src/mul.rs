use super::Gf2Poly;
use super::Limb;
use crate::LimbStorage;
use crate::BITS;

unsafe extern "C" {
    fn gf2x_mul(
        c: *mut Limb,
        a: *const Limb,
        an: Limb,
        b: *const Limb,
        bn: Limb,
    ) -> core::ffi::c_int;
}

impl core::ops::Mul<&Gf2Poly> for &Gf2Poly {
    type Output = Gf2Poly;

    fn mul(self, rhs: &Gf2Poly) -> Self::Output {
        if self.is_zero() || rhs.is_zero() {
            return Gf2Poly::default();
        }

        let capacity = self.limbs().len() + rhs.limbs().len();
        let mut limbs = LimbStorage::with_capacity(capacity);
        let ret = unsafe {
            // safety: gf2x requires that the return buffer has an + bn
            // capacity.
            gf2x_mul(
                limbs.as_mut_ptr(),
                self.limbs().as_ptr(),
                self.limbs().len() as Limb,
                rhs.limbs().as_ptr(),
                rhs.limbs().len() as Limb,
            )
        };
        match ret {
            0 => {}
            -2 => panic!("Out of memory while multiplying gf2 polynomials"),
            _ => panic!("Error while multiplying gf2 polynomials"),
        }

        let deg = self.deg() + rhs.deg();
        unsafe {
            limbs.set_len((deg / BITS + 1) as usize);
        }

        Gf2Poly { limbs, deg }
    }
}

impl core::ops::MulAssign<&Gf2Poly> for Gf2Poly {
    fn mul_assign(&mut self, rhs: &Gf2Poly) {
        *self = &*self * rhs;
    }
}

impl core::ops::Mul<Gf2Poly> for Gf2Poly {
    type Output = Self;

    fn mul(self, rhs: Gf2Poly) -> Self::Output {
        &self * &rhs
    }
}

impl core::ops::Mul<Gf2Poly> for &Gf2Poly {
    type Output = Gf2Poly;

    fn mul(self, rhs: Gf2Poly) -> Self::Output {
        self * &rhs
    }
}

impl core::ops::Mul<&Gf2Poly> for Gf2Poly {
    type Output = Gf2Poly;

    fn mul(self, rhs: &Gf2Poly) -> Self::Output {
        &self * rhs
    }
}

#[cfg(test)]
mod tests {
    use crate::prop_assert_poly_eq;

    use super::*;
    use proptest::prelude::*;
    #[test]
    fn mul() {
        let a: Gf2Poly = "5".parse().unwrap();
        let b: Gf2Poly = "6".parse().unwrap();
        let c = a * b;
        assert_eq!(c.to_string(), "1e");
    }

    #[test]
    fn mul_big() {
        let a: Gf2Poly = "f85ad5ccb96b8764b39e85da047ff33a".parse().unwrap();
        let b: Gf2Poly = "156fd09d3fb6f2e4abb2bcbd80a0f47e".parse().unwrap();
        let c = a * b;
        assert_eq!(
            c.to_string(),
            "cb8e7aa6c53847c815aa4bae5e67613776eb938a2e82940a5b257ecec89a12c"
        );
    }

    proptest! {
        #[test]
        fn associativity(a: Gf2Poly, b: Gf2Poly, c: Gf2Poly) {
            prop_assert_poly_eq!(&(&a * &b) * &c, a * (b * c));
        }

        #[test]
        fn commutativity(a: Gf2Poly, b: Gf2Poly) {
            prop_assert_poly_eq!(&a * &b, b * a);
        }

        #[test]
        fn identity(a: Gf2Poly) {
            prop_assert_poly_eq!(&a * &Gf2Poly::one(), a);
        }

        #[test]
        fn distributivity(a: Gf2Poly, b: Gf2Poly, c: Gf2Poly) {
            prop_assert_poly_eq!(&a * &(&b + &c), &a * &b + &a * &c);
        }

        #[test]
        fn degree_homo(a: Gf2Poly, b: Gf2Poly) {
            prop_assume!(!a.is_zero() && !b.is_zero());
            prop_assert_eq!(a.deg() + b.deg(), (a * b).deg());
        }
    }
}
