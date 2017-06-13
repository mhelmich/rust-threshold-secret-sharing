// Copyright (c) 2016 rust-threshold-secret-sharing developers
//
// Licensed under the Apache License, Version 2.0
// <LICENSE-APACHE or http://www.apache.org/licenses/LICENSE-2.0> or the MIT
// license <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. All files in the project carrying such notice may not be copied,
// modified, or distributed except according to those terms.

use rand;
use std::borrow::Borrow;

/// Abstract (Finite) Field definition.
///
/// This trait is not meant to represent a general field in the strict
/// mathematical sense but it has everything we need to make the algorithms work.
pub trait Field 
{
    type E;
    
    fn zero(&self) -> Self::E;
    
    fn one(&self) -> Self::E;
        
    fn add<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E;
    
    fn sub<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E;
    
    fn mul<A: Borrow<Self::E>, B: Borrow<Self::E>>(&self, a: A, b: B) -> Self::E;
    
    fn pow<A: Borrow<Self::E>>(&self, a: A, e: u32) -> Self::E;
    
    fn inv<A: Borrow<Self::E>>(&self, a: A) -> Self::E;
    
    fn eq<L: Borrow<Self::E>, R: Borrow<Self::E>>(&self, lhs: L, rhs: R) -> bool;
    
    fn neq<L: Borrow<Self::E>, R: Borrow<Self::E>>(&self, lhs: L, rhs: R) -> bool {
        ! self.eq(lhs, rhs)
    }
    
    fn sample_with_replacement<R: rand::Rng>(&self, count: usize, rng: &mut R) -> Vec<Self::E>;
}

pub trait PrimeField : Field
{
    type P;
}

pub trait New<T>
where Self: Field
{
    fn new(params: T) -> Self;
}

pub trait Encode<T>
where Self: Field
{
    fn encode(&self, x: T) -> Self::E;
}

pub trait Decode<U>
where Self: Field
{
    fn decode<E: Borrow<Self::E>>(&self, e: E) -> U;
}

/// Helper trait for encoding values to field elements.
pub trait SliceEncode<T>
where Self: Field
{
    fn encode_slice<V: AsRef<[T]>>(&self, values: V) -> Vec<Self::E>;
}

impl<F, T> SliceEncode<T> for F 
where F: Field + Encode<T>, T: Copy
{
    fn encode_slice<V: AsRef<[T]>>(&self, values: V) -> Vec<Self::E> {
        values.as_ref().iter().map(|&v| self.encode(v)).collect()
    }
}

/// Helper trait for decoding field elements to values.
pub trait SliceDecode<U>
where Self: Field
{
    fn decode_slice<E: AsRef<[Self::E]>>(&self, elements: E) -> Vec<U>;
}

impl<F, U> SliceDecode<U> for F
where F: Field + Decode<U>, F::E: Clone
{
    fn decode_slice<E: AsRef<[F::E]>>(&self, elements: E) -> Vec<U> {
        elements.as_ref().iter().map(|e| self.decode(e)).collect()
    }
}

#[allow(unused_macros)]
macro_rules! all_fields_test {
    ($field:ty) => {
        #[test] fn test_convert() { ::fields::test::test_convert::<$field>(); }
        #[test] fn test_add() { ::fields::test::test_add::<$field>(); }
        #[test] fn test_sub() { ::fields::test::test_sub::<$field>(); }
        #[test] fn test_mul() { ::fields::test::test_mul::<$field>(); }
        #[test] fn test_pow() { ::fields::test::test_pow::<$field>(); }
        #[test] fn test_fft2() { ::numtheory::fft::test::test_fft2::<$field>(); }
        #[test] fn test_fft2_inverse() { ::numtheory::fft::test::test_fft2_inverse::<$field>(); }
        #[test] fn test_fft2_big() { ::numtheory::fft::test::test_fft2_big::<$field>(); }
        #[test] fn test_fft3() { ::numtheory::fft::test::test_fft3::<$field>(); }
        #[test] fn test_fft3_inverse() { ::numtheory::fft::test::test_fft3_inverse::<$field>(); }
        #[test] fn test_fft3_big() { ::numtheory::fft::test::test_fft3_big::<$field>(); }
    }
}

#[cfg(test)]
pub mod test {
    use super::*;
    
    pub fn test_convert<F>() 
    where F: Field + PrimeField + New<u32> + Encode<u32> + Decode<u32>, F::P: From<u32>
    {
        let zp = F::new(17);
        for i in 0_u32..20 {
            assert_eq!(zp.decode(zp.encode(i)), i % 17);
        }
    }

    pub fn test_add<F>() 
    where F: Field + PrimeField + New<u32> + Encode<u32> + Decode<u32>, F::P: From<u32>
    {
        let zp = F::new(17);
        assert_eq!(zp.decode(zp.add(zp.encode(8), zp.encode(2))), 10);
        assert_eq!(zp.decode(zp.add(zp.encode(8), zp.encode(13))), 4);
    }

    pub fn test_sub<F>() 
    where F: Field + PrimeField + New<u32> + Encode<u32> + Decode<u32>, F::P: From<u32>
    {
        let zp = F::new(17);
        assert_eq!(zp.decode(zp.sub(zp.encode(8), zp.encode(2))), 6);
        assert_eq!(zp.decode(zp.sub(zp.encode(8), zp.encode(13))), 12);
    }

    pub fn test_mul<F>() 
    where F: Field + PrimeField + New<u32> + Encode<u32> + Decode<u32>, F::P: From<u32>
    {
        let zp = F::new(17);
        assert_eq!(zp.decode(zp.mul(zp.encode(8), zp.encode(2))), (8 * 2) % 17);
        assert_eq!(zp.decode(zp.mul(zp.encode(8), zp.encode(5))), (8 * 5) % 17);
    }

    pub fn test_pow<F>() 
    where F: Field + PrimeField + New<u32> + Encode<u32> + Decode<u32>, F::P: From<u32>
    {
        let zp = F::new(17);
        assert_eq!(zp.decode(zp.pow(zp.encode(2), 0)), 1);
        assert_eq!(zp.decode(zp.pow(zp.encode(2), 3)), 8);
        assert_eq!(zp.decode(zp.pow(zp.encode(2), 6)), 13);
    }
}


mod natural;
pub use self::natural::NaturalPrimeField;

mod montgomery;
pub use self::montgomery::MontgomeryField32;

#[cfg(feature="largefield")] mod large;
#[cfg(feature="largefield")] pub use self::large::LargePrimeField;

// pub mod native;
