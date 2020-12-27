//! Testing tools which are not bound to a particular module

use std::panic::{catch_unwind, UnwindSafe};

/// Check that a function panics when called
pub fn panics<R>(f: impl FnOnce() -> R + UnwindSafe) -> bool {
    catch_unwind(f).is_err()
}

/// Check that a certain user input passes maximally rigorous validation
///
/// Here, we use the concept of a user input validation function, which is a
/// recuring pattern in this crate. This function can either panic if the input
/// is fatally wrong, or return false if it is a bit weird. Sometimes, unit
/// tests need to take a given action whether an input is weird or wrong.
///
pub fn is_standard<I: UnwindSafe>(
    input: I,
    validation: impl FnOnce(I) -> bool + UnwindSafe,
) -> bool {
    match catch_unwind(|| validation(input)) {
        Ok(true) => true,
        Ok(false) | Err(_) => false,
    }
}
