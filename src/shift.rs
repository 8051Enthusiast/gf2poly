use core::ops::Bound;
use core::ops::RangeBounds;

use crate::BITS;
use crate::LimbStorage;
use crate::limbs_for_deg;

use crate::Gf2Poly;

macro_rules! shift_impl {
    ($op:tt, $trait:ident, $fun:ident, $impl:ty, $($t:ty),+) => {
        $(
            impl core::ops::$trait<$t> for $impl {
                type Output = Gf2Poly;

                fn $fun(self, rhs: $t) -> Self::Output {
                    self $op rhs as u64
                }
            }
        )+
    };
}

macro_rules! shift_assign_impl {
    ($op:tt, $trait:ident, $fun:ident, $($t:ty),+) => {
        $(
            impl core::ops::$trait<$t> for Gf2Poly {
                fn $fun(&mut self, rhs: $t) {
                    *self $op rhs as u64;
                }
            }
        )+
    };
}

impl core::ops::ShlAssign<u64> for Gf2Poly {
    fn shl_assign(&mut self, rhs: u64) {
        if rhs == 0 || self.is_zero() {
            return;
        }
        let bit_shift = rhs % BITS;
        let Ok(limb_shift) = usize::try_from(rhs / BITS) else {
            panic!("Left shift of {rhs} is too big")
        };
        let new_deg = self
            .deg
            .checked_add(rhs)
            .expect("Left shift would overflow degree");
        let old_size = self.limbs().len();
        self.resize_to_deg_unnormalized(new_deg);

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
    }
}

impl core::ops::Shl<u64> for Gf2Poly {
    type Output = Self;

    fn shl(mut self, rhs: u64) -> Self::Output {
        self <<= rhs;
        self
    }
}

shift_impl!(<<, Shl, shl, Gf2Poly, usize, u32, u16, u8);
shift_assign_impl!(<<=, ShlAssign, shl_assign, usize, u32, u16, u8);

impl core::ops::ShrAssign<u64> for Gf2Poly {
    fn shr_assign(&mut self, rhs: u64) {
        if rhs == 0 || self.is_zero() {
            return;
        }

        let bit_shift = rhs % BITS;
        let Ok(limb_shift) = usize::try_from(rhs / BITS) else {
            *self = Self::zero();
            return;
        };
        let Some(new_deg) = self.deg.checked_sub(rhs) else {
            *self = Self::zero();
            return;
        };
        let limb_size = self.limbs().len();

        if bit_shift == 0 {
            self.limbs.copy_within(limb_shift.., 0);
            self.limbs.truncate(limb_size - limb_shift);
            self.deg = new_deg;
        } else {
            for i in limb_shift..limb_size {
                let src = self.limbs[i];
                let src_val = src;
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
    }
}

impl core::ops::Shr<u64> for Gf2Poly {
    type Output = Self;

    fn shr(mut self, rhs: u64) -> Self::Output {
        self >>= rhs;
        self
    }
}

shift_impl!(>>, Shr, shr, Gf2Poly, usize, u32, u16, u8);
shift_assign_impl!(>>=, ShrAssign, shr_assign, usize, u32, u16, u8);

impl Gf2Poly {
    fn normalize_range<T: RangeBounds<u64>>(&self, range: T) -> Option<(u64, u64)> {
        if self.is_zero() {
            return None;
        }

        let lo = match range.start_bound() {
            Bound::Included(n) => *n,
            Bound::Excluded(n) => n.checked_add(1)?,
            Bound::Unbounded => 0,
        };

        let hi = match range.end_bound() {
            Bound::Included(n) => *n,
            Bound::Excluded(n) => n.checked_sub(1)?,
            Bound::Unbounded => self.deg(),
        };

        let hi = self.deg().min(hi);

        if lo > hi {
            return None;
        }

        Some((lo, hi))
    }

    /// Creates a copy of the bit range given in `range`.
    ///
    /// ## Example
    /// ```rust
    /// # use gf2poly::Gf2Poly;
    /// let p: Gf2Poly = "123456789".parse().unwrap();
    /// let subrange = p.subrange(12..25);
    /// assert_eq!(subrange.to_string(), "1456");
    /// ```
    pub fn subrange<T: RangeBounds<u64>>(&self, range: T) -> Gf2Poly {
        let Some((lo, hi)) = self.normalize_range(range) else {
            return Gf2Poly::zero();
        };

        let bit_shift = lo % BITS;
        let limb_start = (lo / BITS) as usize;
        let limb_end = (hi / BITS + 1) as usize;
        let limb_count = limbs_for_deg(hi - lo);

        let mut limbs;
        if bit_shift == 0 {
            limbs = LimbStorage::from(&self.limbs()[limb_start..limb_end]);
        } else {
            limbs = LimbStorage::with_capacity(limb_count);
            for i in limb_start..limb_end {
                let src = self.limbs[i];
                if let Some(target) = limbs.get_mut((i - limb_start).wrapping_sub(1)) {
                    *target |= src << (BITS - bit_shift);
                };
                if limbs.len() < limb_count {
                    limbs.push(src >> bit_shift);
                }
            }
        }

        let highest_bit_offset = (hi - lo) % BITS;
        *limbs.last_mut().unwrap() &= (1 << highest_bit_offset) | ((1 << highest_bit_offset) - 1);
        Gf2Poly::from_limb_storage(limbs)
    }

    /// Equivalent to *self += rhs << shift, but more efficient
    /// since rhs << shift is not calculated explicitely.
    pub fn fused_shl_add(&mut self, rhs: &Gf2Poly, shift: u64) {
        if rhs.is_zero() {
            return;
        }
        let bit_shift = shift % BITS;
        let Ok(limb_shift) = usize::try_from(shift / BITS) else {
            panic!("Left shift of {rhs} is too big")
        };
        let shift_deg = rhs
            .deg
            .checked_add(shift)
            .expect("Left shift would overflow degree");
        let new_deg = shift_deg.max(self.deg());
        self.resize_to_deg_unnormalized(new_deg);

        if bit_shift == 0 {
            for (self_limb, rhs_limb) in self.limbs[limb_shift..].iter_mut().zip(rhs.limbs.iter()) {
                *self_limb ^= *rhs_limb;
            }
        } else {
            let rhs_count = rhs.limbs.len();
            for i in (0..rhs_count).rev() {
                let src = rhs.limbs[i];
                if let Some(target) = self.limbs.get_mut(limb_shift + i + 1) {
                    *target ^= src >> (BITS - bit_shift);
                };
                self.limbs[limb_shift + i] ^= src << bit_shift;
            }
            if matches!(self.limbs().last(), Some(0)) {
                self.limbs.pop();
            }
        }
        self.normalize();
    }
}

impl core::ops::Shr<u64> for &Gf2Poly {
    type Output = Gf2Poly;

    fn shr(self, rhs: u64) -> Self::Output {
        self.subrange(rhs..)
    }
}

shift_impl!(>>, Shr, shr, &Gf2Poly, usize, u32, u16, u8);

impl core::ops::Shl<u64> for &Gf2Poly {
    type Output = Gf2Poly;

    fn shl(self, rhs: u64) -> Self::Output {
        let copy: Gf2Poly = self.clone();
        copy << rhs
    }
}

shift_impl!(<<, Shl, shl, &Gf2Poly, usize, u32, u16, u8);

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

        #[test]
        fn shift_impl_same(a: Gf2Poly, n in 0u32..256) {
            prop_assert_poly_eq!(a.clone() >> n, &a >> n);
            prop_assert_poly_eq!(a.clone() << n, &a << n);
        }

        #[test]
        fn split_and_recombine(a: Gf2Poly, idx1 in 0..128u64, idx2 in 0..128u64) {
            prop_assume!(!a.is_constant());
            let mut idx1 = idx1 % a.deg();
            let mut idx2 = idx2 % a.deg();
            if idx2 < idx1 {
                core::mem::swap(&mut idx1, &mut idx2);
            }

            let lo = a.subrange(..idx1);
            let mid = a.subrange(idx1..idx2);
            let hi = a.subrange(idx2..=a.deg());

            let recombined = lo + mid * Gf2Poly::x_to_the_power_of(idx1) + hi * Gf2Poly::x_to_the_power_of(idx2);

            prop_assert_poly_eq!(a, recombined);
        }

        #[test]
        fn fused_shl_add(a: Gf2Poly, b: Gf2Poly, shift in 0u64..256) {
            let mut copy = a.clone();
            copy.fused_shl_add(&b, shift);
            let shifted = a + (b << shift);
            prop_assert_poly_eq!(copy, shifted);
        }
    }
}
