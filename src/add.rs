// lol characteristic 2
#![allow(clippy::suspicious_arithmetic_impl)]
use crate::{LimbStorage, deg, normalized_len};

use super::Gf2Poly;

impl core::ops::Add<&Gf2Poly> for &Gf2Poly {
    type Output = Gf2Poly;

    fn add(self, rhs: &Gf2Poly) -> Self::Output {
        let new_size = self.limbs().len().max(rhs.limbs().len());
        let mut limbs = LimbStorage::with_capacity(new_size);
        limbs.extend(
            self.limbs()
                .iter()
                .zip(rhs.limbs().iter())
                .map(|(a, b)| a ^ b),
        );
        if rhs.limbs().len() == self.limbs().len() {
            // only in this case can the highest terms cancel out
            limbs.truncate(normalized_len(&limbs));
        } else if limbs.len() < self.limbs().len() {
            // the zip only wrote the common limbs, so here we add the rest
            // from whichever polynomial was bigger
            limbs.extend_from_slice(&self.limbs()[limbs.len()..]);
        } else {
            limbs.extend_from_slice(&rhs.limbs()[limbs.len()..]);
        }
        let deg = deg(&limbs);
        Gf2Poly { limbs, deg }
    }
}

impl core::ops::AddAssign<&Gf2Poly> for Gf2Poly {
    fn add_assign(&mut self, rhs: &Gf2Poly) {
        let new_size = self.limbs().len().max(rhs.limbs().len());
        self.limbs.resize(new_size, 0);
        for (a, b) in self.limbs.iter_mut().zip(rhs.limbs().iter()) {
            *a ^= b;
        }
        if self.deg == rhs.deg {
            self.limbs.truncate(normalized_len(&self.limbs))
        } else if rhs.deg() > self.deg() {
            // in the case self.deg() > rhs.deg(), self already has the right
            // upper bits, so we only have to handle this case
            self.limbs
                .extend_from_slice(&rhs.limbs()[self.limbs.len()..]);
        }
        self.deg = deg(&self.limbs);
    }
}

impl core::ops::AddAssign<Gf2Poly> for Gf2Poly {
    fn add_assign(&mut self, mut rhs: Gf2Poly) {
        if rhs.limbs().len() > self.limbs().len() {
            core::mem::swap(self, &mut rhs);
        }
        *self += &rhs;
    }
}

impl core::ops::Add<Gf2Poly> for Gf2Poly {
    type Output = Gf2Poly;

    fn add(mut self, mut rhs: Gf2Poly) -> Self::Output {
        if self.limbs().len() < rhs.limbs().len() {
            core::mem::swap(&mut self, &mut rhs);
        }
        self += &rhs;
        self
    }
}

impl core::ops::Add<Gf2Poly> for &Gf2Poly {
    type Output = Gf2Poly;

    fn add(self, rhs: Gf2Poly) -> Self::Output {
        self.clone() + rhs
    }
}

impl core::ops::Add<&Gf2Poly> for Gf2Poly {
    type Output = Gf2Poly;

    fn add(self, rhs: &Gf2Poly) -> Self::Output {
        self + rhs.clone()
    }
}

impl core::ops::Neg for Gf2Poly {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self
    }
}

impl core::ops::Neg for &Gf2Poly {
    type Output = Gf2Poly;

    fn neg(self) -> Self::Output {
        self.clone()
    }
}

impl core::ops::Sub<Gf2Poly> for Gf2Poly {
    type Output = Gf2Poly;

    fn sub(self, rhs: Gf2Poly) -> Self::Output {
        self + rhs
    }
}

impl core::ops::Sub<Gf2Poly> for &Gf2Poly {
    type Output = Gf2Poly;

    fn sub(self, rhs: Gf2Poly) -> Self::Output {
        self + rhs
    }
}

impl core::ops::Sub<&Gf2Poly> for Gf2Poly {
    type Output = Gf2Poly;

    fn sub(self, rhs: &Gf2Poly) -> Self::Output {
        self + rhs
    }
}

impl core::ops::Sub<&Gf2Poly> for &Gf2Poly {
    type Output = Gf2Poly;

    fn sub(self, rhs: &Gf2Poly) -> Self::Output {
        self + rhs
    }
}
#[cfg(test)]
mod tests {
    use crate::prop_assert_poly_eq;

    use super::*;
    use proptest::prelude::*;
    #[test]
    fn add() {
        let a: Gf2Poly = "5".parse().unwrap();
        let b: Gf2Poly = "6".parse().unwrap();
        let c = a + b;
        assert_eq!(c.to_string(), "3");
    }

    proptest! {
        #[test]
        fn associativity(a: Gf2Poly, b: Gf2Poly, c: Gf2Poly) {
            prop_assert_poly_eq!(&(&a + &b) + &c, a + (b + c));
        }

        #[test]
        fn commutativity(a: Gf2Poly, b: Gf2Poly) {
            prop_assert_poly_eq!(&a + &b, b + a);
        }

        #[test]
        fn identity(a: Gf2Poly) {
            let zero = Gf2Poly::zero();
            prop_assert_poly_eq!(&a + &zero, a);
        }

        #[test]
        fn inverse(a: Gf2Poly) {
            prop_assert_poly_eq!(&a + -&a, Gf2Poly::zero());
        }

        #[test]
        fn char2(a: Gf2Poly) {
            prop_assert_poly_eq!(&a + &a, Gf2Poly::zero());
        }

        #[test]
        fn add_assign(mut a: Gf2Poly, b: Gf2Poly) {
            let c = &a + &b;
            a += &b;
            prop_assert_poly_eq!(a, c);
        }

        #[test]
        fn add_assign_consume(mut a: Gf2Poly, b: Gf2Poly) {
            let c = &a + &b;
            a += b;
            prop_assert_poly_eq!(a, c);
        }
    }
}
