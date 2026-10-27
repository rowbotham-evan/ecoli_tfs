// Loads data/dataset.csv and renders a searchable, sortable table.
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

let header = [], rows = [], sortCol = -1, sortDir = 1, methodCol = -1, methodFilter = "";

// Hidden and compressed columns are remembered per browser, by column name.
const VIEW_KEY = "column-view";
const columnView = loadView();

function loadView() {
  try {
    const v = JSON.parse(localStorage.getItem(VIEW_KEY)) || {};
    return { hidden: new Set(v.hidden || []), compressed: new Set(v.compressed || []) };
  } catch { return { hidden: new Set(), compressed: new Set() }; }
}

function saveView() {
  try {
    localStorage.setItem(VIEW_KEY, JSON.stringify({ hidden: [...columnView.hidden], compressed: [...columnView.compressed] }));
  } catch {}
}

function renderColumnToggles() {
  document.getElementById("columns").innerHTML = "Columns: " + header.map(h =>
    `<label><input type="checkbox" data-col="${escape(h)}"${columnView.hidden.has(h) ? "" : " checked"}> ${escape(h)}</label>`).join("");
}

// The method column's header is a dropdown that filters rows by method.
function headerCell(h, i) {
  const arrow = i === sortCol ? (sortDir > 0 ? " ▲" : " ▼") : "";
  const compressed = columnView.compressed.has(h);
  const toggle = `<button type="button" class="compress" data-col="${escape(h)}" title="${compressed ? "Expand" : "Compress"} column">${compressed ? "+" : "−"}</button>`;
  const attrs = `data-i="${i}" title="${escape(h)}"${compressed ? ' class="compressed"' : ""}`;
  if (i !== methodCol || compressed) return `<th ${attrs}>${toggle}${escape(h)}${arrow}</th>`;
  const options = [["", `${h}: all`], ...METHODS.map(m => [m, m])];
  return `<th ${attrs}>${toggle}<select id="method-filter" aria-label="Filter by method">` + options.map(([v, label]) =>
    `<option value="${v}"${v === methodFilter ? " selected" : ""}>${escape(label)}</option>`).join("") + `</select>${arrow}</th>`;
}

function bodyCell(v, i) {
  const h = header[i];
  if (columnView.compressed.has(h)) return `<td class="compressed" title="${escape(v)}">${escape(v)}</td>`;
  if (h.endsWith("_sequence")) return `<td class="seq">${escape(v)}</td>`;
  return `<td title="${escape(v)}">${escape(v.length > 60 ? v.slice(0, 60) + "…" : v)}</td>`;
}

function render() {
  const q = search.value.trim().toLowerCase();
  let shown = rows.filter(r =>
    (!methodFilter || r[methodCol] === methodFilter) &&
    (!q || r.some(v => v.toLowerCase().includes(q))));
  if (sortCol >= 0) {
    shown.sort((a, b) => sortDir * a[sortCol].localeCompare(b[sortCol], undefined, { numeric: true }));
  }
  const cols = header.map((h, i) => i).filter(i => !columnView.hidden.has(header[i]));
  table.tHead.innerHTML = "<tr>" + cols.map(i => headerCell(header[i], i)).join("") + "</tr>";
  const body = shown.slice(0, MAX_ROWS).map(r =>
    "<tr>" + cols.map(i => bodyCell(r[i], i)).join("") + "</tr>");
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
    return;
  }
  const th = e.target.closest("th");
  if (!th || e.target.closest("select")) return;
  const i = Number(th.dataset.i);
  sortDir = i === sortCol ? -sortDir : 1;
  sortCol = i;
  render();
});
table.tHead.addEventListener("change", e => {
  if (e.target.id !== "method-filter") return;
  methodFilter = e.target.value;
  render();
});
search.addEventListener("input", render);
document.getElementById("columns").addEventListener("change", e => {
  const h = e.target.dataset.col;
  if (!h) return;
  e.target.checked ? columnView.hidden.delete(h) : columnView.hidden.add(h);
  saveView();
  render();
});

fetch(DATA_URL)
  .then(r => { if (!r.ok) throw new Error(r.status); return r.text(); })
  .then(text => {
    [header = [], ...rows] = parseCSV(text);
    if (!header.length) throw new Error("empty");
    methodCol = header.indexOf("method");
    renderColumnToggles();
    render();
    empty.hidden = rows.length > 0;
  })
  .catch(() => {
    empty.hidden = false;
    table.hidden = true;
    search.disabled = true;
  });
