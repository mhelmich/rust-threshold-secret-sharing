// Copyright (c) 2017 rust-threshold-secret-sharing developers

//! Prime field using GMP's Mpz as the underlying type.

extern crate gmp;

use std::borrow::Borrow;
use rand;

use ::fields::{Field, PrimeField, New, Encode, Decode};
use ::numtheory::{generic_mod_pow};

#[derive(Clone,Debug,PartialEq)]
pub struct LargePrimeField(gmp::mpz::Mpz);

impl Field for LargePrimeField
{
    /// Invariant is that numbers are stored in canonical form [0..prime).
    type E = gmp::mpz::Mpz;

    fn zero(&self) -> Self::E {
        gmp::mpz::Mpz::zero()
    }

    fn one(&self) -> Self::E {
        gmp::mpz::Mpz::one()
    }

    fn add<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E {
        (a.borrow() + b.borrow()) % &self.0
    }

    fn sub<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E {
        let c = (a.borrow() - b.borrow()) % &self.0;
        if c >= gmp::mpz::Mpz::zero() {
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
        let (mut s, mut old_s) = (gmp::mpz::Mpz::zero(), gmp::mpz::Mpz::one());
        let (mut t, mut old_t) = (gmp::mpz::Mpz::one(),  gmp::mpz::Mpz::zero());
        let (mut r, mut old_r) = (self.0.clone(), a.borrow().clone());

        let mut tmp = gmp::mpz::Mpz::zero();
        while ! r.is_zero() {
            let quotient = &old_r / &r;
            tmp.clone_from(&r); r = &old_r - &quotient * r; old_r.clone_from(&tmp);
            tmp.clone_from(&s); s = &old_s - &quotient * s; old_s.clone_from(&tmp);
            tmp.clone_from(&t); t = &old_t - &quotient * t; old_t.clone_from(&tmp);
        }

        let d = old_r;
        let inv = old_s % &self.0;
        debug_assert_eq!(d, gmp::mpz::Mpz::one());
        if inv >= gmp::mpz::Mpz::zero() {
            inv
        } else {
            inv + &self.0
        }
    }
    
    fn eq<L: Borrow<Self::E>, R: Borrow<Self::E>>(&self, lhs: L, rhs: R) -> bool {
        lhs.borrow() == rhs.borrow()
    }
    
    fn sample_with_replacement<R: rand::Rng>(&self, count: usize, rng: &mut R) -> Vec<Self::E> {
        // TODO verify and comment
        let upper = &self.0;
        let bitsize = upper.bit_length();
        let bytes = (bitsize-1) / 8 + 1;
        (0..count).map(|_| {
            loop {
                let mut buf: Vec<u8> = vec![0; bytes];
                rng.fill_bytes(&mut buf);
                let n = gmp::mpz::Mpz::from(&*buf) >> (bytes*8-bitsize);
                if n < *upper {
                    return n
                }
            }
        }).collect()
    }
}

impl PrimeField for LargePrimeField {
    type P = gmp::mpz::Mpz;
}


impl New<gmp::mpz::Mpz> for LargePrimeField {
    fn new(prime: gmp::mpz::Mpz) -> Self {
        LargePrimeField(prime)
    }
}

impl<'a> New<&'a gmp::mpz::Mpz> for LargePrimeField {
    fn new(prime: &'a gmp::mpz::Mpz) -> Self {
        Self::new(prime.clone())
    }
}

impl<'a> New<&'a str> for LargePrimeField {    
    fn new(prime: &'a str) -> Self {
        use std::str::FromStr;
        Self::new(gmp::mpz::Mpz::from_str(&prime).unwrap())
    }
}

impl New<u32> for LargePrimeField 
{
    fn new(prime: u32) -> Self {
        Self::new(gmp::mpz::Mpz::from(prime))
    }
}

impl New<u64> for LargePrimeField 
{
    fn new(prime: u64) -> Self {
        Self::new(gmp::mpz::Mpz::from(prime))
    }
}


impl<'a> Encode<&'a gmp::mpz::Mpz> for LargePrimeField
{
    fn encode(&self, x: &'a gmp::mpz::Mpz) -> Self::E {
        let y = x % &self.0;
        if y >= gmp::mpz::Mpz::zero() {
            y
        } else {
            y + &self.0
        }
    }
}

impl Encode<gmp::mpz::Mpz> for LargePrimeField
{
    fn encode(&self, x: gmp::mpz::Mpz) -> Self::E {
        self.encode(&x)
    }
}

impl<'a> Encode<&'a str> for LargePrimeField
{    
    fn encode(&self, x: &'a str) -> Self::E {
        use std::str::FromStr;
        self.encode(gmp::mpz::Mpz::from_str(&x).unwrap())
    }
}

impl Encode<u32> for LargePrimeField 
{
    fn encode(&self, x: u32) -> Self::E {
        self.encode(gmp::mpz::Mpz::from(x))
    }
}

impl Encode<u64> for LargePrimeField 
{
    fn encode(&self, x: u64) -> Self::E {
        self.encode(gmp::mpz::Mpz::from(x))
    }
}

impl Decode<u64> for LargePrimeField
{
    fn decode<E: Borrow<Self::E>>(&self, x: E) -> u64 {
        let foo: Option<u64> = x.borrow().into();
        foo.unwrap() // TODO
    }
}

impl Decode<u32> for LargePrimeField
{
    fn decode<E: Borrow<Self::E>>(&self, x: E) -> u32 {
        let foo: u64 = self.decode(x);
        foo as u32
    }
}

#[cfg(test)]
all_fields_test!(LargePrimeField);
