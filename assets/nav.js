// Navigation bar links open as soon as the mouse button goes down, without waiting for
// it to be released. Only a plain left click is handled: clicks with a modifier key
// (new tab, new window), other buttons, the keyboard and touch keep their normal behavior.
document.addEventListener("mousedown", e => {
  if (e.button !== 0 || e.metaKey || e.ctrlKey || e.shiftKey || e.altKey) return;
  const link = e.target.closest(".site-title a, .site-tabs a");
  if (!link) return;
  e.preventDefault();
  location.href = link.href;
});

// Swimming E. coli along the whole bar (styled and animated in assets/nav.css).
// Each cell is a small SVG: a rod with three flagella. A flagellum is a long wavy line
// that slides along itself behind a clipping window, which reads as a wave travelling
// down it. Purely decorative, so it is hidden from screen readers.
(() => {
  const bar = document.querySelector(".site-bar-inner");
  if (!bar) return;
  // fill, outline, size, vertical position, seconds per crossing (of the old middle strip), seconds between tumbles
  const CELLS = [
    ["#9fe0bd", "#3f9a6e", 1.00, 30, 26, 11],
    ["#a9d8ee", "#3d8db3", 0.88, 68, 35, 14],
    ["#fbd3a4", "#d98632", 0.80, 26, 21, 9],
    ["#cfc5f3", "#7b68c9", 0.94, 72, 30, 13],
    ["#f6bfd0", "#c9587f", 0.72, 48, 41, 17],
  ];
  const wave = "M-24 22" + " q4 -5 8 0 t8 0".repeat(6);   // 16 units per wavelength, longer than its window
  const strip = document.createElement("div");
  strip.className = "site-bugs";
  strip.setAttribute("aria-hidden", "true");
  const LONGER = 2.4;   // the full bar is about this many times the old middle strip: keep the same swimming speed
  strip.innerHTML = CELLS.map(([fill, line, size, top, crossing, tumble], i) => {
    const run = crossing * LONGER;
    const flagellum = angle => `<g transform="rotate(${angle} 44 22)" clip-path="url(#site-bug-clip-${i})"><path class="wave" d="${wave}"/></g>`;
    // The negative delays start every cell part-way through its run; "left" is where it rests without animation.
    return `<div class="site-bug" style="top: ${top}%; width: ${Math.round(96 * size)}px; left: ${8 + i * 19}%; animation-duration: ${run}s; animation-delay: ${-(run * (0.13 + i * 0.21) % run).toFixed(1)}s">` +
      `<div class="site-bug-bob" style="animation-delay: ${-i * 0.7}s"><div class="site-bug-tumble" style="animation-duration: ${tumble}s; animation-delay: ${-i * 2.3}s">` +
      `<svg viewBox="0 0 96 44"><clipPath id="site-bug-clip-${i}"><rect x="2" y="6" width="42" height="32"/></clipPath>` +
      `<g fill="none" stroke="${line}" stroke-width="1.5" stroke-linecap="round" opacity="0.8">${flagellum(0)}${flagellum(17)}${flagellum(-17)}</g>` +
      `<rect x="42" y="12" width="48" height="20" rx="10" fill="${fill}" stroke="${line}" stroke-width="1.6"/>` +
      `<path d="M54 22 q5 -6 10 0 t10 0" fill="none" stroke="${line}" stroke-width="1.4" stroke-linecap="round" opacity="0.55"/>` +
      `<circle cx="51" cy="19" r="1.3" fill="${line}" opacity="0.5"/><circle cx="81" cy="26" r="1.3" fill="${line}" opacity="0.5"/>` +
      `</svg></div></div></div>`;
  }).join("");
  bar.appendChild(strip);

  // The cells swim the whole bar and turn almost transparent behind its labels. The mask is
  // a left-to-right gradient: opaque in the gaps, faint over the title, each tab and the
  // download button, with a short fade at every edge. Measured again whenever the bar's
  // layout can change (resize, fonts arriving).
  const FAINT = "rgba(0, 0, 0, 0.12)", FADE = 14, PAD = 8;
  const setMask = () => {
    const box = strip.getBoundingClientRect();
    if (!box.width) return;
    const spans = [...bar.querySelectorAll(".site-title, .site-tabs a, .site-download")]
      .map(el => el.getBoundingClientRect()).filter(r => r.width)
      .map(r => [Math.max(0, r.left - box.left - PAD), Math.min(box.width, r.right - box.left + PAD)])
      .sort((a, b) => a[0] - b[0]);
    const merged = [];   // labels too close for a cell to show between them count as one
    for (const span of spans) {
      const last = merged[merged.length - 1];
      if (last && span[0] - last[1] < 2 * FADE + 6) last[1] = Math.max(last[1], span[1]);
      else merged.push(span);
    }
    const px = v => `${Math.round(Math.min(Math.max(v, 0), box.width))}px`;
    const stops = merged.flatMap(([a, b]) => [`#000 ${px(a - FADE)}`, `${FAINT} ${px(a)}`, `${FAINT} ${px(b)}`, `#000 ${px(b + FADE)}`]);
    strip.style.setProperty("--bugs-mask", `linear-gradient(90deg, ${["#000 0px", ...stops, "#000 100%"].join(", ")})`);
  };
  setMask();
  window.addEventListener("resize", setMask);
  window.addEventListener("load", setMask);
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(setMask);
})();
