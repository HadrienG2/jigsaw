//! Avoid duplicate logger initialization in tests

use std::sync::Once;

/// Initialize env_logger once
pub fn init_logger() {
    static LOGGER_INIT: Once = Once::new();
    LOGGER_INIT.call_once(|| env_logger::init());
}
