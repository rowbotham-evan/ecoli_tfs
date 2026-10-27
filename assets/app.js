// Loads data/dataset.csv and renders a searchable table.
const DATA_URL = "data/dataset.csv";
const MAX_ROWS = 500;
const METHODS = ["SELEX", "TTPT", "DAP", "CHIP"];

function parseCSV(text) {
  const rows = [];
  let row = [], field = "", quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') { field += '"'; i++; }
      else if (c === '"') quoted = false;
      else field += c;
    } else if (c === '"') quoted = true;
    else if (c === ",") { row.push(field); field = ""; }
    else if (c === "\n" || c === "\r") {
      if (c === "\r" && text[i + 1] === "\n") i++;
      row.push(field); rows.push(row); row = []; field = "";
    } else field += c;
  }
  if (field || row.length) { row.push(field); rows.push(row); }
  return rows.filter(r => r.some(v => v !== ""));
}

const table = document.getElementById("table");
const search = document.getElementById("search");
const count = document.getElementById("count");
const empty = document.getElementById("empty");

let header = [], rows = [], methodCol = -1, methodFilter = "";

// Compressed columns are remembered per browser, by column name.
const VIEW_KEY = "column-view";
const columnView = loadView();

function loadView() {
  try {
    const v = JSON.parse(localStorage.getItem(VIEW_KEY)) || {};
    return { compressed: new Set(v.compressed || []) };
  } catch { return { compressed: new Set() }; }
}

function saveView() {
  try {
    localStorage.setItem(VIEW_KEY, JSON.stringify({ compressed: [...columnView.compressed] }));
  } catch {}
}

function headerCell(h, i) {
  const compressed = columnView.compressed.has(h);
  const toggle = `<button type="button" class="compress" data-col="${escape(h)}" title="${compressed ? "Expand" : "Compress"} column">${compressed ? "+" : "−"}</button>`;
  const attrs = `data-i="${i}" title="${escape(h)}"${compressed ? ' class="compressed"' : ""}`;
  if (i !== methodCol || compressed) return `<th ${attrs}>${toggle}${escape(h)}</th>`;
  const options = [["", `${h}: all`], ...METHODS.map(m => [m, m])];
  return `<th ${attrs}>${toggle}<select id="method-filter" aria-label="Filter by method">` + options.map(([v, label]) =>
    `<option value="${v}"${v === methodFilter ? " selected" : ""}>${escape(label)}</option>`).join("") + `</select></th>`;
}

function bodyCell(v, i) {
  const h = header[i];
  if (columnView.compressed.has(h)) return `<td class="compressed" title="${escape(v)}">${escape(v)}</td>`;
  if (h.endsWith("_sequence")) return `<td class="seq">${escape(v)}</td>`;
  return `<td title="${escape(v)}">${escape(v.length > 60 ? v.slice(0, 60) + "…" : v)}</td>`;
}

function render() {
  const q = search.value.trim().toLowerCase();
  const shown = rows.filter(r =>
    (!methodFilter || r[methodCol] === methodFilter) &&
    (!q || r.some(v => v.toLowerCase().includes(q))));
  table.tHead.innerHTML = "<tr>" + header.map(headerCell).join("") + "</tr>";
  const body = shown.slice(0, MAX_ROWS).map(r =>
    "<tr>" + r.map(bodyCell).join("") + "</tr>");
  table.tBodies[0].innerHTML = body.join("");
  count.textContent = `${shown.length} of ${rows.length} entries` + (shown.length > MAX_ROWS ? ` (showing first ${MAX_ROWS})` : "");
}

function escape(s) {
  return s.replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

table.tHead.addEventListener("click", e => {
  const button = e.target.closest("button.compress");
  if (button) {
    const h = button.dataset.col;
    columnView.compressed.has(h) ? columnView.compressed.delete(h) : columnView.compressed.add(h);
    saveView();
    render();
  }
});
table.tHead.addEventListener("change", e => {
  if (e.target.id !== "method-filter") return;
  methodFilter = e.target.value;
  render();
});
search.addEventListener("input", render);

fetch(DATA_URL)
  .then(r => { if (!r.ok) throw new Error(r.status); return r.text(); })
  .then(text => {
    [header = [], ...rows] = parseCSV(text);
    if (!header.length) throw new Error("empty");
    methodCol = header.indexOf("method");
    render();
    empty.hidden = rows.length > 0;
  })
  .catch(() => {
    empty.hidden = false;
    table.hidden = true;
    search.disabled = true;
  });
