// Dev-only auto-reload: polls watched files and reloads the page when one changes.
// Only runs on localhost, so it does nothing on GitHub Pages.
(() => {
  if (!["localhost", "127.0.0.1"].includes(location.hostname)) return;
  const files = window.WATCH_FILES || [location.pathname];
  const stamps = new Map();
  let first = true;

  async function check() {
    for (const f of files) {
      try {
        const r = await fetch(f, { method: "HEAD", cache: "no-store" });
        const t = r.headers.get("Last-Modified") + "|" + r.headers.get("Content-Length");
        if (!first && stamps.get(f) !== t) { location.reload(); return; }
        stamps.set(f, t);
      } catch {}
    }
    first = false;
    setTimeout(check, 1000);
  }
  check();
})();
