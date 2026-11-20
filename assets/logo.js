// Motif logo drawing shared by index.html (the TF dashboard) and motif.html.
// Exposes window.Logo. The *Svg functions return a bare <svg>; the block functions
// wrap it in a <figure> with an enlarge button, a caption and a table of the numbers.
window.Logo = (() => {
  const BASES = ["A", "C", "G", "T"];
  // Conventional nucleotide colors (A green, C blue, G orange, T red).
  const BASE_COLORS = { A: "#1f9e5a", C: "#2b6cb0", G: "#e09b12", T: "#d23b3b" };
  const LOGO_FONT = "700 100px Helvetica, Arial, sans-serif";
  const esc = s => String(s).replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const glyphBox = (() => {
    const ctx = document.createElement("canvas").getContext("2d");
    ctx.font = LOGO_FONT;
    return Object.fromEntries(BASES.map(b => {
      const m = ctx.measureText(b);
      return [b, { left: m.actualBoundingBoxLeft, w: m.actualBoundingBoxLeft + m.actualBoundingBoxRight,
                   desc: m.actualBoundingBoxDescent, h: m.actualBoundingBoxAscent + m.actualBoundingBoxDescent }];
    }));
  })();
  // One letter stretched to fill a box colW wide and h tall whose bottom edge is at yBottom.
  const glyph = (b, x, yBottom, h, colW) => {
    const g = glyphBox[b], gx = x + colW * 0.06, gw = colW * 0.88;
    return `<text transform="translate(${gx.toFixed(2)} ${yBottom.toFixed(2)}) scale(${(gw / g.w).toFixed(4)} ${(h / g.h).toFixed(4)}) translate(${g.left.toFixed(2)} ${(-g.desc).toFixed(2)})" ` +
      `style="font: ${LOGO_FONT}" fill="${BASE_COLORS[b]}">${b}</text>`;
  };

  // Per-position base frequencies of equal-length aligned sequences.
  function frequencies(seqs) {
    return [...seqs[0]].map((_, i) => {
      const f = { A: 0, C: 0, G: 0, T: 0 };
      for (const s of seqs) f[s[i]] += 1 / seqs.length;
      return f;
    });
  }
  // Rows of [A, C, G, T] numbers (counts or probabilities) as per-position frequencies.
  function fromMatrix(matrix) {
    return matrix.map(row => {
      const total = row.reduce((a, b) => a + b, 0) || 1;
      return Object.fromEntries(BASES.map((b, i) => [b, row[i] / total]));
    });
  }
  const bits = f => 2 + BASES.reduce((sum, b) => sum + (f[b] > 0 ? f[b] * Math.log2(f[b]) : 0), 0);

  // A logo the reader can click to see enlarged in the page's #zoom dialog. The caption
  // is HTML built by the page; it is stored escaped in data-caption and shown as HTML.
  function zoomButton(media, title, captionHtml) {
    return `<button type="button" class="zoom" data-title="${esc(title)}" data-caption="${esc(captionHtml)}" aria-label="Enlarge: ${esc(title)}">` +
      `${media}<span class="zoom-hint" aria-hidden="true">Enlarge</span></button>`;
  }

  // Sequence logo: letter height = base frequency x information content (bits) at
  // that position, largest letter on top.
  function sequenceSvg(freqs) {
    const colW = 30, left = 40, top = 10, plotH = 120, bottom = 22;
    const W = left + freqs.length * colW + 6, H = top + plotH + bottom;
    const y = v => top + plotH - (v / 2) * plotH;
    const parts = [];
    for (const t of [0, 1, 2]) {
      parts.push(`<line x1="${left}" x2="${W - 6}" y1="${y(t)}" y2="${y(t)}" stroke="${t ? "#ecebe6" : "#c3c2b7"}" stroke-width="1"/>`);
      parts.push(`<text x="${left - 8}" y="${y(t)}" text-anchor="end" dominant-baseline="central" font-size="11" fill="#666">${t}</text>`);
    }
    parts.push(`<text transform="translate(12 ${top + plotH / 2}) rotate(-90)" text-anchor="middle" font-size="11" fill="#666">bits</text>`);
    freqs.forEach((f, i) => {
      const ic = bits(f);
      const x = left + i * colW;
      let yBottom = y(0);
      for (const b of [...BASES].sort((p, q) => f[p] - f[q])) {      // smallest at the bottom
        const h = f[b] * ic / 2 * plotH;
        if (h < 0.5) continue;
        parts.push(glyph(b, x, yBottom, h, colW));
        yBottom -= h;
      }
      parts.push(`<text x="${x + colW / 2}" y="${H - 6}" text-anchor="middle" font-size="10" fill="#666">${i + 1}</text>`);
      const tip = `Position ${i + 1}: ` + BASES.map(b => `${b} ${f[b].toFixed(2)}`).join(" · ") + ` · ${ic.toFixed(2)} bits`;
      parts.push(`<rect x="${x}" y="${top}" width="${colW}" height="${plotH}" fill="transparent"><title>${tip}</title></rect>`);
    });
    return `<svg viewBox="0 0 ${W} ${H}" style="max-width: ${W}px" role="img" aria-label="Sequence logo, ${freqs.length} positions">${parts.join("")}</svg>`;
  }

  function sequence(freqs, caption, title) {
    const table = `<details><summary>View as table</summary><table><thead><tr><th>Position</th>${BASES.map(b => `<th>${b}</th>`).join("")}<th>Bits</th></tr></thead><tbody>` +
      freqs.map((f, i) => `<tr><td>${i + 1}</td>${BASES.map(b => `<td>${f[b].toFixed(2)}</td>`).join("")}<td>${bits(f).toFixed(2)}</td></tr>`).join("") +
      `</tbody></table></details>`;
    const fullCaption = `${caption} Letter height: frequency × information content.`;
    return `<figure class="logo">${zoomButton(sequenceSvg(freqs), title, fullCaption)}<figcaption>${fullCaption}</figcaption>${table}</figure>`;
  }

  // Energy (weight) logo: per position, positive weights stack upward from zero and
  // negative weights downward, largest furthest from zero; letter height = |weight|.
  function energySvg(matrix) {
    const colW = 18, left = 60, top = 10, plotH = 150, bottom = 22;   // about the published logos' proportions
    const L = matrix.length, W = left + L * colW + 6, H = top + plotH + bottom;
    const up = Math.max(0.1, ...matrix.map(row => row.reduce((s, v) => s + Math.max(v, 0), 0)));
    const dn = Math.max(0.1, ...matrix.map(row => row.reduce((s, v) => s + Math.max(-v, 0), 0)));
    const scale = plotH / (up + dn), y0 = top + up * scale;
    const step = up + dn > 2 ? 0.5 : 0.25;
    const parts = [];
    for (let t = -Math.floor(dn / step) * step; t <= up + 1e-9; t += step) {
      const yy = y0 - t * scale;
      parts.push(`<line x1="${left}" x2="${W - 6}" y1="${yy.toFixed(2)}" y2="${yy.toFixed(2)}" stroke="${Math.abs(t) < 1e-9 ? "#c3c2b7" : "#ecebe6"}" stroke-width="1"/>`);
      parts.push(`<text x="${left - 8}" y="${yy.toFixed(2)}" text-anchor="end" dominant-baseline="central" font-size="11" fill="#666">${(Math.abs(t) < 1e-9 ? 0 : t).toFixed(step < 0.5 ? 2 : 1)}</text>`);
    }
    parts.push(`<text transform="translate(12 ${top + plotH / 2}) rotate(-90)" text-anchor="middle" font-size="11" fill="#666">weight</text>`);
    matrix.forEach((row, i) => {
      const x = left + i * colW;
      const cells = BASES.map((b, k) => [b, row[k]]);
      let yb = y0;
      for (const [b, v] of cells.filter(c => c[1] > 0).sort((p, q) => p[1] - q[1])) {
        const h = v * scale;
        if (h >= 0.5) parts.push(glyph(b, x, yb, h, colW));
        yb -= h;
      }
      let yt = y0;
      for (const [b, v] of cells.filter(c => c[1] < 0).sort((p, q) => q[1] - p[1])) {
        const h = -v * scale;
        if (h >= 0.5) parts.push(glyph(b, x, yt + h, h, colW));
        yt += h;
      }
      parts.push(`<text x="${x + colW / 2}" y="${H - 6}" text-anchor="middle" font-size="10" fill="#666">${i + 1}</text>`);
      const tip = `Position ${i + 1}: ` + cells.map(([b, v]) => `${b} ${v >= 0 ? "+" : ""}${v.toFixed(2)}`).join(" · ");
      parts.push(`<rect x="${x}" y="${top}" width="${colW}" height="${plotH}" fill="transparent"><title>${tip}</title></rect>`);
    });
    return `<svg viewBox="0 0 ${W} ${H}" style="max-width: ${W}px" role="img" aria-label="Energy logo, ${L} positions">${parts.join("")}</svg>`;
  }

  function energy(matrix, bias, caption, title) {
    const table = `<details><summary>View as table</summary><table><thead><tr><th>Position</th>${BASES.map(b => `<th>${b}</th>`).join("")}</tr></thead><tbody>` +
      matrix.map((row, i) => `<tr><td>${i + 1}</td>${row.map(v => `<td>${v.toFixed(3)}</td>`).join("")}</tr>`).join("") +
      `</tbody></table><p class="zoom-cap">Bias term: ${bias.toFixed(3)}</p></details>`;
    const fullCaption = `${caption} Above zero favors binding; letter height is the weight.`;
    return `<figure class="logo">${zoomButton(energySvg(matrix), title, fullCaption)}<figcaption>${fullCaption}</figcaption>${table}</figure>`;
  }

  return { BASES, frequencies, fromMatrix, bits, zoomButton, sequenceSvg, sequence, energySvg, energy };
})();
