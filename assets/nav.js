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

// Swimming E. coli in the bar's empty middle (styled and animated in assets/nav.css).
// Each cell is a small SVG: a rod with three flagella. A flagellum is a long wavy line
// that slides along itself behind a clipping window, which reads as a wave travelling
// down it. Purely decorative, so it is hidden from screen readers.
(() => {
  const bar = document.querySelector(".site-bar-inner");
  if (!bar) return;
  // fill, outline, size, vertical position, seconds per crossing, seconds between tumbles
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
  strip.innerHTML = CELLS.map(([fill, line, size, top, run, tumble], i) => {
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
})();
