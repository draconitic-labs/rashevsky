//! Desktop workbench of the Rashevsky engine. The window shows `desktop/ui`; this process holds the
//! open models and runs the analyses on background threads.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod commands;
mod jobs;
mod payload;

fn main() {
    // This machine's WebKitGTK crashes on Wayland with its DMA-BUF renderer ("Error 71"), and the
    // same flags keep software canvases smooth elsewhere. Anything set by the caller wins.
    #[cfg(target_os = "linux")]
    for (key, value) in [("WEBKIT_DISABLE_DMABUF_RENDERER", "1"), ("GDK_BACKEND", "x11"), ("WEBKIT_DISABLE_COMPOSITING_MODE", "1")] {
        if std::env::var_os(key).is_none() {
            // SAFETY: nothing else runs yet, so no thread can read the environment concurrently.
            unsafe { std::env::set_var(key, value) };
        }
    }
    let _ = rayon::ThreadPoolBuilder::new().stack_size(64 << 20).build_global();
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .manage(commands::AppState::default())
        .invoke_handler(tauri::generate_handler![
            commands::startup,
            commands::open_bundled,
            commands::open_file,
            commands::download_model,
            commands::solve,
            commands::reaction_range,
            commands::layout,
            commands::start_screen,
            commands::start_sampling,
            commands::poll_job,
            commands::cancel_job,
            commands::envelope,
            commands::variability,
            commands::kinetic_open,
            commands::kinetic_simulate,
            commands::kinetic_control,
            commands::kinetic_stochastic,
            commands::atlas,
            commands::export_fluxes,
            commands::export_escher,
            commands::export_sbml,
            commands::validation,
            commands::qa_log,
        ])
        .run(tauri::generate_context!())
        .expect("the window failed to open");
}
