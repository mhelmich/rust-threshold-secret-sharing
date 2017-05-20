// Copyright (c) 2017 rust-threshold-secret-sharing developers

#[macro_use]
extern crate bencher;
extern crate ramp;

use bencher::Bencher;

pub fn bench_add_ref_ref(bencher: &mut Bencher) 
{    
    bencher.iter(|| {
        let ref a = ramp::Int::from(1231231);
        let ref b = ramp::Int::from(423421);
        let _ = a + b;
    });
}

pub fn bench_add_ref_val(bencher: &mut Bencher) 
{        
    bencher.iter(|| {
        let ref a = ramp::Int::from(1231231);
        let     b = ramp::Int::from(423421);
        let _ = a + b;
    });
}

pub fn bench_add_val_ref(bencher: &mut Bencher) 
{    
    bencher.iter(|| {
        let     a = ramp::Int::from(1231231);
        let ref b = ramp::Int::from(423421);
        let _ = a + b;
    });
}

pub fn bench_add_val_val(bencher: &mut Bencher) 
{    
    bencher.iter(|| {
        let a = ramp::Int::from(1231231);
        let b = ramp::Int::from(423421);
        let _ = a + b;
    });
}

benchmark_group!(add,
    bench_add_ref_ref,
    bench_add_ref_val,
    bench_add_val_ref,
    bench_add_val_val);

benchmark_main!(add);
