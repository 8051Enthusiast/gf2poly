use crate::BITS;

use crate::Gf2Poly;

impl core::ops::Shl<u64> for Gf2Poly {
    type Output = Self;

    fn shl(mut self, rhs: u64) -> Self::Output {
        if rhs == 0 {
            return self;
        }
        if self.is_zero() {
            return self;
        }
        let Ok(shift) = usize::try_from(rhs) else {
            panic!("Left shift of {rhs} is too big")
        };
        let new_deg = self
            .deg
            .checked_add(rhs)
            .expect("Left shift would overflow degree");
        let old_size = self.limbs().len();
        self.resize_to_deg_unnormalized(new_deg);

        let limb_shift = shift / BITS;
        let bit_shift = shift % BITS;
        if bit_shift == 0 {
            self.limbs.copy_within(..old_size, limb_shift);
            self.limbs[..limb_shift].fill(0);
        } else {
            for i in (0..old_size).rev() {
                let src = &mut self.limbs[i];
                let src_val = *src;
                *src = 0;
                if let Some(target) = self.limbs.get_mut(limb_shift + i + 1) {
                    *target |= src_val >> (BITS - bit_shift);
                };
                self.limbs[limb_shift + i] = src_val << bit_shift;
            }
            if matches!(self.limbs().last(), Some(0)) {
                self.limbs.pop();
            }
        }
        self
    }
}

impl core::ops::Shl<u32> for Gf2Poly {
    type Output = Self;

    fn shl(self, rhs: u32) -> Self::Output {
        self << rhs as u64
    }
}

impl core::ops::Shl<usize> for Gf2Poly {
    type Output = Self;

    fn shl(self, rhs: usize) -> Self::Output {
        self << rhs as u64
    }
}

impl core::ops::Shl<u16> for Gf2Poly {
    type Output = Self;

    fn shl(self, rhs: u16) -> Self::Output {
        self << rhs as u64
    }
}

impl core::ops::Shl<u8> for Gf2Poly {
    type Output = Self;

    fn shl(self, rhs: u8) -> Self::Output {
        self << rhs as u64
    }
}

impl core::ops::ShlAssign<u64> for Gf2Poly {
    fn shl_assign(&mut self, rhs: u64) {
        *self = core::mem::take(self) << rhs;
    }
}

impl core::ops::ShlAssign<u32> for Gf2Poly {
    fn shl_assign(&mut self, rhs: u32) {
        *self = core::mem::take(self) << rhs;
    }
}

impl core::ops::ShlAssign<usize> for Gf2Poly {
    fn shl_assign(&mut self, rhs: usize) {
        *self = core::mem::take(self) << rhs;
    }
}

impl core::ops::ShlAssign<u16> for Gf2Poly {
    fn shl_assign(&mut self, rhs: u16) {
        *self = core::mem::take(self) << rhs;
    }
}

impl core::ops::ShlAssign<u8> for Gf2Poly {
    fn shl_assign(&mut self, rhs: u8) {
        *self = core::mem::take(self) << rhs;
    }
}

impl core::ops::Shr<u64> for Gf2Poly {
    type Output = Self;

    fn shr(mut self, rhs: u64) -> Self::Output {
        if rhs == 0 {
            return self;
        }
        if self.is_zero() {
            return self;
        }
        let Ok(shift) = usize::try_from(rhs) else {
            return Self::zero();
        };
        let Some(new_deg) = self.deg.checked_sub(rhs) else {
            return Self::zero();
        };
        let limb_size = self.limbs().len();

        let limb_shift = shift / BITS;
        let bit_shift = shift % BITS;
        if bit_shift == 0 {
            self.limbs.copy_within(limb_shift.., 0);
            self.limbs.truncate(limb_size - limb_shift);
            self.deg = new_deg;
        } else {
            for i in limb_shift..limb_size {
                let src = &mut self.limbs[i];
                let src_val = *src;
                if let Some(target) = self.limbs.get_mut((i - limb_shift).wrapping_sub(1)) {
                    *target |= src_val << (BITS - bit_shift);
                };
                self.limbs[i - limb_shift] = src_val >> bit_shift;
            }
            self.resize_to_deg_unnormalized(new_deg);
            if matches!(self.limbs.last(), Some(0)) {
                self.limbs.pop();
            }
        }
        self
    }
}

impl core::ops::Shr<u32> for Gf2Poly {
    type Output = Self;

    fn shr(self, rhs: u32) -> Self::Output {
        self >> rhs as u64
    }
}

impl core::ops::Shr<usize> for Gf2Poly {
    type Output = Self;

    fn shr(self, rhs: usize) -> Self::Output {
        self >> rhs as u64
    }
}

impl core::ops::Shr<u16> for Gf2Poly {
    type Output = Self;

    fn shr(self, rhs: u16) -> Self::Output {
        self >> rhs as u64
    }
}

impl core::ops::Shr<u8> for Gf2Poly {
    type Output = Self;

    fn shr(self, rhs: u8) -> Self::Output {
        self >> rhs as u64
    }
}

impl core::ops::ShrAssign<u64> for Gf2Poly {
    fn shr_assign(&mut self, rhs: u64) {
        *self = core::mem::take(self) >> rhs;
    }
}

impl core::ops::ShrAssign<u32> for Gf2Poly {
    fn shr_assign(&mut self, rhs: u32) {
        *self = core::mem::take(self) >> rhs;
    }
}

impl core::ops::ShrAssign<usize> for Gf2Poly {
    fn shr_assign(&mut self, rhs: usize) {
        *self = core::mem::take(self) >> rhs;
    }
}

impl core::ops::ShrAssign<u16> for Gf2Poly {
    fn shr_assign(&mut self, rhs: u16) {
        *self = core::mem::take(self) >> rhs;
    }
}

impl core::ops::ShrAssign<u8> for Gf2Poly {
    fn shr_assign(&mut self, rhs: u8) {
        *self = core::mem::take(self) >> rhs;
    }
}

#[cfg(test)]
mod tests {
    use crate::prop_assert_poly_eq;

    use super::*;
    use proptest::prelude::*;
    #[test]
    fn shift_left() {
        let a: Gf2Poly = "dbf".parse().unwrap();
        let b = a.clone() << 3u32;
        assert_eq!(b.to_string(), "6df8");
        let b = a.clone() << 51u32;
        assert_eq!(b.to_string(), "6df8000000000000");
        let b = a.clone() << 63u32;
        assert_eq!(b.to_string(), "6df8000000000000000");
        let b = a.clone() << 64u32;
        assert_eq!(b.to_string(), "dbf0000000000000000");
        let b = a.clone() << 65u32;
        assert_eq!(b.to_string(), "1b7e0000000000000000");
    }

    #[test]
    fn shift_right() {
        let a: Gf2Poly = "dbf39fa3902865daff3".parse().unwrap();
        let b = a.clone() >> 3u32;
        assert_eq!(b.to_string(), "1b7e73f472050cbb5fe");
        let b = a.clone() >> 51u32;
        assert_eq!(b.to_string(), "1b7e73f");
        let b = a.clone() >> 63u32;
        assert_eq!(b.to_string(), "1b7e");
        let b = a.clone() >> 64u32;
        assert_eq!(b.to_string(), "dbf");
        let b = a.clone() >> 65u32;
        assert_eq!(b.to_string(), "6df");
    }

    proptest! {
        #[test]
        fn shift_parts_left(a: Gf2Poly, n in 0u32..256, m in 0u32..256) {
            prop_assert_poly_eq!(((a.clone() << n) << m) >> (n + m), a);
        }

        #[test]
        fn shift_parts_right(a: Gf2Poly, n in 0u32..256, m in 0u32..256) {
            prop_assert_poly_eq!(((a.clone() << (n + m)) >> n) >> m, a);
        }

        #[test]
        fn shift_right_deg(a: Gf2Poly) {
            prop_assume!(!a.is_zero());
            let deg = a.deg();
            prop_assert_poly_eq!(a >> deg, Gf2Poly::one());
        }
    }
}
