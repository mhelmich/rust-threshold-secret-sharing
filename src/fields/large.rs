// Copyright (c) 2017 rust-threshold-secret-sharing developers

//! Prime field using RAMP as the underlying type.

extern crate ramp;

use std::borrow::Borrow;
use rand;
use rand::Rng;
use self::ramp::RandomInt;

use ::fields::Field;
use ::fields::PrimeField;
use ::fields::{Encode, Decode};
use ::numtheory::{generic_mod_pow};

#[derive(Clone,Debug,PartialEq)]
pub struct LargePrimeField(ramp::Int);

// impl LargePrimeField {
//     fn egcd(a: &Self, b: &Self) -> (Self, Self, Self) {
//         if b == 0 {
//             (a.clone(), Self::one(), Self::zero())
//         } else {
//             let q = a / b;
//             let r = a % b;
//             let (d, s, t) = Self::egcd(b, &r);
//             let new_t = s - &t * q;
//             (d, t, new_t)
//         }
//     }
// }

impl Field for LargePrimeField
{
    /// Invariant is that numbers are stored in canonical form in [0..prime).
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
        let (mut t, mut old_t) = (ramp::Int::one(),  ramp::Int::zero());
        let (mut r, mut old_r) = (self.0.clone(), a.borrow().clone());

        let mut tmp = ramp::Int::zero();
        while r != 0 {
            let quotient = &old_r / &r;
            tmp.clone_from(&r); r = &old_r - &quotient * r; old_r.clone_from(&tmp);
            tmp.clone_from(&s); s = &old_s - &quotient * s; old_s.clone_from(&tmp);
            tmp.clone_from(&t); t = &old_t - &quotient * t; old_t.clone_from(&tmp);
        }

        let d = old_r;
        let inv = old_s % &self.0;
        assert_eq!(d, ramp::Int::one());
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
        (0..count).map(|_| rng.gen_uint_below(&self.0)).collect()
        // unimplemented!();
    }
}

impl PrimeField for LargePrimeField {
    type P = ramp::Int;
    
    fn new(prime: Self::P) -> Self {
        LargePrimeField(prime)
    }    
}

impl<T> Encode<T> for LargePrimeField 
where ramp::Int: From<T>
{
    fn encode(&self, x: T) -> Self::E {
        let y = ramp::Int::from(x) % &self.0;
        if y >= 0 {
            y
        } else {
            y + &self.0
        }
    }
}

impl Decode<u32> for LargePrimeField {
    fn decode<E: Borrow<Self::E>>(&self, x: E) -> u32 {
        u32::from(x.borrow())
    }
}

#[cfg(test)]
all_fields_test!(LargePrimeField);
