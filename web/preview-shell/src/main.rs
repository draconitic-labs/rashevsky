//! Serves `web/dist` in a Tauri window, so the demo is tested in the webview it targets.
fn main() {
    for (key, value) in
        [("WEBKIT_DISABLE_DMABUF_RENDERER", "1"), ("GDK_BACKEND", "x11"), ("WEBKIT_DISABLE_COMPOSITING_MODE", "1")]
    {
        if std::env::var_os(key).is_none() {
            // SAFETY: nothing else runs yet.
            unsafe { std::env::set_var(key, value) };
        }
    }
    tauri::Builder::default().run(tauri::generate_context!()).expect("the window failed to open");
}
