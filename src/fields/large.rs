// Copyright (c) 2017 rust-threshold-secret-sharing developers

//! Prime field using RAMP as the underlying type.

extern crate framp as ramp;

use rand;
use std::borrow::Borrow;

use fields::{Decode, Encode, Field, New, PrimeField};
use numtheory::generic_mod_pow;

#[derive(Clone, Debug, PartialEq)]
pub struct LargePrimeField(ramp::Int);

impl Field for LargePrimeField {
    /// Invariant is that numbers are stored in canonical form [0..prime).
    type E = ramp::Int;

    fn zero(&self) -> Self::E {
        ramp::Int::zero()
    }

    fn one(&self) -> Self::E {
        ramp::Int::one()
    }

    fn add<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E {
        (a.borrow() + b.borrow()) % &self.0
    }

    fn sub<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E {
        let c = (a.borrow() - b.borrow()) % &self.0;
        if c >= 0 {
            c
        } else {
            c + &self.0
        }
    }

    fn mul<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E {
        (a.borrow() * b.borrow()) % &self.0
    }

    fn pow<A: Borrow<Self::E>>(&self, a: A, e: u32) -> Self::E {
        generic_mod_pow(self, a.borrow().clone(), e)
    }

    fn inv<A: Borrow<Self::E>>(&self, a: A) -> Self::E {
        let (mut s, mut old_s) = (ramp::Int::zero(), ramp::Int::one());
        let (mut t, mut old_t) = (ramp::Int::one(), ramp::Int::zero());
        let (mut r, mut old_r) = (self.0.clone(), a.borrow().clone());

        let mut tmp = ramp::Int::zero();
        while r != 0 {
            let quotient = &old_r / &r;
            tmp.clone_from(&r);
            r = &old_r - &quotient * r;
            old_r.clone_from(&tmp);
            tmp.clone_from(&s);
            s = &old_s - &quotient * s;
            old_s.clone_from(&tmp);
            tmp.clone_from(&t);
            t = &old_t - &quotient * t;
            old_t.clone_from(&tmp);
        }

        let d = old_r;
        let inv = old_s % &self.0;
        debug_assert_eq!(d, ramp::Int::one());
        if inv >= 0 {
            inv
        } else {
            inv + &self.0
        }
    }

    fn eq<L: Borrow<Self::E>, R: Borrow<Self::E>>(&self, lhs: L, rhs: R) -> bool {
        lhs.borrow() == rhs.borrow()
    }

    fn sample_with_replacement<R: rand::Rng>(&self, count: usize, rng: &mut R) -> Vec<Self::E> {
        use self::ramp::RandomInt;
        (0..count).map(|_| rng.gen_uint_below(&self.0)).collect()
    }
}

impl PrimeField for LargePrimeField {
    type P = ramp::Int;
}

impl New<ramp::Int> for LargePrimeField {
    fn new(prime: ramp::Int) -> Self {
        LargePrimeField(prime)
    }
}

impl<'a> New<&'a ramp::Int> for LargePrimeField {
    fn new(prime: &'a ramp::Int) -> Self {
        Self::new(prime.clone())
    }
}

impl<'a> New<&'a str> for LargePrimeField {
    fn new(prime: &'a str) -> Self {
        use std::str::FromStr;
        Self::new(ramp::Int::from_str(&prime).unwrap())
    }
}

impl New<usize> for LargePrimeField {
    fn new(prime: usize) -> Self {
        Self::new(ramp::Int::from(prime))
    }
}

impl New<u8> for LargePrimeField {
    fn new(prime: u8) -> Self {
        Self::new(ramp::Int::from(prime))
    }
}

impl New<u16> for LargePrimeField {
    fn new(prime: u16) -> Self {
        Self::new(ramp::Int::from(prime))
    }
}

impl New<u32> for LargePrimeField {
    fn new(prime: u32) -> Self {
        Self::new(ramp::Int::from(prime))
    }
}

impl New<u64> for LargePrimeField {
    fn new(prime: u64) -> Self {
        Self::new(ramp::Int::from(prime))
    }
}

impl<'a> Encode<&'a ramp::Int> for LargePrimeField {
    fn encode(&self, x: &'a ramp::Int) -> Self::E {
        let y = x % &self.0;
        if y >= 0 {
            y
        } else {
            y + &self.0
        }
    }
}

impl Encode<ramp::Int> for LargePrimeField {
    fn encode(&self, x: ramp::Int) -> Self::E {
        self.encode(&x)
    }
}

impl<'a> Encode<&'a str> for LargePrimeField {
    fn encode(&self, x: &'a str) -> Self::E {
        use std::str::FromStr;
        self.encode(ramp::Int::from_str(&x).unwrap())
    }
}

impl Encode<usize> for LargePrimeField {
    fn encode(&self, x: usize) -> Self::E {
        self.encode(ramp::Int::from(x))
    }
}

impl Encode<u8> for LargePrimeField {
    fn encode(&self, x: u8) -> Self::E {
        self.encode(ramp::Int::from(x))
    }
}

impl Encode<u16> for LargePrimeField {
    fn encode(&self, x: u16) -> Self::E {
        self.encode(ramp::Int::from(x))
    }
}

impl Encode<u32> for LargePrimeField {
    fn encode(&self, x: u32) -> Self::E {
        self.encode(ramp::Int::from(x))
    }
}

impl Encode<u64> for LargePrimeField {
    fn encode(&self, x: u64) -> Self::E {
        self.encode(ramp::Int::from(x))
    }
}

// impl<T> Encode<T> for LargePrimeField
// where ramp::Int: From<T>
// {
//     fn encode(&self, x: T) -> Self::E {
//         self.encode(ramp::Int::from(x))
//     }
// }

impl<U> Decode<U> for LargePrimeField
where
    for<'a> U: From<&'a ramp::Int>,
{
    fn decode<E: Borrow<Self::E>>(&self, x: E) -> U {
        U::from(x.borrow())
    }
}

#[cfg(test)]
all_fields_test!(LargePrimeField);
