// How well a TF's two published motifs agree. Shared by coverage.html (all TFs) and
// motif.html (one TF). Exposes window.MotifAgreement.
//
// The scores are not computed here: scripts/build_motif_agreement.py aligns each pair of
// motifs with Tomtom (MEME Suite), allowing shifts and reverse complements, and scores
// the aligned position probability matrices (PPMs) by Pearson correlation. It writes
// data/motif_agreement.json, which this file loads.
//
// load() resolves to { tfs, mismatched, threshold, method }:
//   tfs         one entry per TF with two motifs, best score first: name, r (Pearson
//               correlation), overlap (aligned positions), reversed (one motif reverse-
//               complemented), tomtom_p (Tomtom's p-value for the alignment), p (share of
//               mismatched pairs scoring as high), other (what the DAP-seq motif was
//               compared with) and points (the aligned cells: x = the base's probability
//               in the DAP-seq motif, y = its probability in the other motif)
//   mismatched  sorted scores of one TF's DAP-seq motif against another TF's other motif
//   threshold   the 95th percentile of mismatched
//   method      the Tomtom version and options, and how energy matrices become PPMs
window.MotifAgreement = (() => {
  const URL = "../data/motif_agreement.json";   // relative to the pages in html/

  function load() {
    return fetch(URL, { cache: "no-store" })
      .then(r => { if (!r.ok) throw new Error(r.status); return r.json(); })
      .then(data => {
        for (const t of data.tfs) t.points = t.points.map(([x, y, k, pos]) => ({ x, y, base: "ACGT"[k], pos }));
        return data;
      });
  }

  // Pearson correlation and least-squares line (y = intercept + slope * x) of a set of points.
  function fit(points) {
    const n = points.length || 1;
    const mx = points.reduce((s, q) => s + q.x, 0) / n, my = points.reduce((s, q) => s + q.y, 0) / n;
    let sxy = 0, sxx = 0, syy = 0;
    for (const q of points) { sxy += (q.x - mx) * (q.y - my); sxx += (q.x - mx) ** 2; syy += (q.y - my) ** 2; }
    const slope = sxx ? sxy / sxx : 0;
    return { r: sxx && syy ? sxy / Math.sqrt(sxx * syy) : 0, slope, intercept: my - slope * mx };
  }

  return { load, fit, URL };
})();
